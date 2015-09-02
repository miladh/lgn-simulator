#include "neuron.h"

Neuron::Neuron(const Config *cfg, Stimuli *stim)
    : m_stim(stim)
{
    const Setting & root = cfg->getRoot();
    m_nPoints = root["spatialDomainSettings"]["nPoints"];

    m_spatialMesh = linspace(-0.5, 0.5, m_nPoints);
    double dr = m_spatialMesh(1) - m_spatialMesh(0);
    double N_2 = ceil(m_nPoints/2.);
    double df = 1./dr/m_nPoints;
    double fs = 1./dr;

    m_freqMesh = linspace(-N_2*df, (m_nPoints - 1. - N_2)*df, m_nPoints);
    m_freqMesh*=2*PI;
//    cout << m_freqMesh << endl;


    m_response = zeros(m_nPoints, m_nPoints);
    m_responseFT = zeros<cx_mat>(m_nPoints, m_nPoints);
    m_impulseResponse = zeros(m_nPoints, m_nPoints);
    m_impulseResponseFT = zeros<cx_mat>(m_nPoints, m_nPoints);


}

Neuron::~Neuron()
{

}


void Neuron::computeResponse(double t)
{
    for(int i = 0; i < m_nPoints; i++){
        for(int j = 0; j < m_nPoints; j++){
            m_responseFT(i,j) = exp(-m_i*m_stim->w() * t)
              *m_stim->frequency({m_freqMesh[i], m_freqMesh[j]}, m_stim->w())
              *impulseResponseFT({m_freqMesh[i], m_freqMesh[j]}, m_stim->w());
        }
    }


    m_responseFT = Functions::fftShift(m_responseFT);
    fftw_complex* in = reinterpret_cast<fftw_complex*> (m_responseFT.memptr());
    fftw_complex* out = reinterpret_cast<fftw_complex*> (m_responseFT.memptr());
    fftw_plan plan2 = fftw_plan_dft_2d(m_nPoints,m_nPoints, in,  out, FFTW_BACKWARD,            FFTW_ESTIMATE);

    fftw_execute(plan2);

    m_responseFT = Functions::fftShift(m_responseFT);
    m_response = real(m_responseFT);

    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////
    ///
//    double dr = m_spatialMesh(1) - m_spatialMesh(0);

//    complexResponse.set_real(m_stim->spatial());
//    complexResponse.set_imag(m_stim->spatial()*0);
//    complexResponse = Functions::fftShift(complexResponse);

//    fftw_complex* in1 = reinterpret_cast<fftw_complex*> (complexResponse.memptr());
//    fftw_complex* out1 = reinterpret_cast<fftw_complex*> (complexResponse.memptr());

//    fftw_plan plan = fftw_plan_dft_2d(complexResponse.n_cols,complexResponse.n_rows,
//                                      in1,  out1, FFTW_FORWARD, FFTW_ESTIMATE);

//    fftw_execute(plan);

//    fftw_plan plan1 = fftw_plan_dft_2d(complexResponse.n_cols,complexResponse.n_rows,
//                                      in1,  out1, FFTW_BACKWARD, FFTW_ESTIMATE);
//    fftw_execute(plan1);

//    complexResponse*=dr*dr;
//    m_stim->setSpatial(Functions::fftShift(real(complexResponse)));

}



void Neuron::computeResponseFT(double w)
{
    cout << "computeResponseComplex: Not implemented!" << endl;
    for(int i = 0; i < m_nPoints; i++){
        for(int j = 0; j < m_nPoints; j++){

            double Gcomplex = impulseResponseFT({m_spatialMesh[i], m_spatialMesh[j]}, w);
            double Scomplex = m_stim->frequency({m_spatialMesh[i], m_spatialMesh[j]}, w);

            m_impulseResponseFT(i,j) = Gcomplex;
            m_responseFT(i,j) = Gcomplex*Scomplex;
        }
    }

}

void Neuron::computeImpulseResponse(double t)
{
    m_impulseResponse  = 0* m_impulseResponse;

//    for(int i = 0; i < int(m_spatialMesh.n_elem); i++){
//        for(int j = 0; j < int(m_spatialMesh.n_elem); j++){
//                        m_impulseResponse(i, j) += 1./8./(PI*PI*PI) *
//                                impulseResponseComplex({x[m], x[n]},x[o])
//                                * w[m] * w[n] * w[o]
//                                * cos(m_spatialMesh[i]*x[m]+ m_spatialMesh[j]*x[n] - x[o]* t);
//                    }

//                }
//            }
}

void Neuron::computeImpulseResponseFT(double w)
{
    cout << "computeImpulseResponseComplex: Not implemented!" << endl;
}



void Neuron::addGanglionCell(Neuron *neuron,
                             SpatialKernel *sKernel,
                             TemporalKernel *tKernel)
{
    m_ganglionCells.emplace_back(Input{neuron, sKernel, tKernel});
}

void Neuron::addRelayCell(Neuron *neuron,
                          SpatialKernel *sKernel,
                          TemporalKernel *tKernel)
{
    m_relayCells.emplace_back(Input{neuron, sKernel, tKernel});
}

void Neuron::addInterNeuron(Neuron *neuron,
                            SpatialKernel *sKernel,
                            TemporalKernel *tKernel)
{
    m_interNeurons.emplace_back(Input{neuron, sKernel, tKernel});
}

void Neuron::addCorticalNeuron(Neuron *neuron,
                               SpatialKernel *sKernel,
                               TemporalKernel *tKernel)
{
    m_corticalNeurons.emplace_back(Input{neuron, sKernel, tKernel});
}


vector<Neuron::Input> Neuron::ganglionCells() const
{
    return m_ganglionCells;
}

vector<Neuron::Input> Neuron::relayCells() const
{
    return m_relayCells;
}


vector<Neuron::Input> Neuron::interNeurons() const
{
    return m_interNeurons;
}

vector<Neuron::Input> Neuron::corticalNeurons() const
{
    return m_corticalNeurons;
}


string Neuron::cellType() const
{
    return m_cellType;
}

mat Neuron::response() const
{
    return m_response;
}

mat Neuron::impulseResponse() const
{
    return m_impulseResponse;
}

cx_mat Neuron::responseFT() const
{
    return m_responseFT;
}

cx_mat Neuron::impulseResponseFT() const
{
    return m_impulseResponseFT;
}










