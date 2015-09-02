#include "neuron.h"

Neuron::Neuron(const Config *cfg, Stimuli *stim)
    : m_stim(stim)
{
    const Setting & root = cfg->getRoot();
    const Setting &gridLimits = root["gridSettings"]["grid"];
    const Setting &integrationDomain = root["gridSettings"]["integrationDomain"];

    int nPoints = root["spatialDomainSettings"]["nPoints"];
//    double dr = root["spatialDomainSettings"]["dr"];


    m_spatialMesh = linspace(-0.5, 0.5, nPoints);
    double dr = m_spatialMesh(1) - m_spatialMesh(0);
    double N_2 = ceil(nPoints/2.);
    double df = 1./dr/nPoints;
    double fs = 1./dr;



    m_freqMesh = linspace(-N_2*df, (nPoints - 1. - N_2)*df, nPoints);
    m_freqMesh*=2*PI;
//    cout << m_freqMesh << endl;

    vec3 grid;
    for(int i = 0; i < 3; i++){
        grid[i] = gridLimits[i];
        m_domain[i] = integrationDomain[i];
    }


    m_responseComplex = zeros(nPoints, nPoints);
    m_impulseResponse = zeros(nPoints, nPoints);
    m_impulseResponseComplex = zeros(nPoints, nPoints);

    complexResponse = zeros<cx_mat>(nPoints, nPoints);
    m_response = zeros(nPoints, nPoints);


}

Neuron::~Neuron()
{

}


void Neuron::computeResponse(double t)
{
    for(int i = 0; i < int(m_freqMesh.n_elem); i++){
        for(int j = 0; j < int(m_freqMesh.n_elem); j++){
            complexResponse(i,j) = exp(-m_i*m_stim->w() * t)
              *m_stim->frequency({m_freqMesh[i], m_freqMesh[j]}, m_stim->w())
              *impulseResponseComplex({m_freqMesh[i], m_freqMesh[j]}, m_stim->w());
        }
    }


    complexResponse = Functions::fftShift(complexResponse);
    fftw_complex* in = reinterpret_cast<fftw_complex*> (complexResponse.memptr());
    fftw_complex* out = reinterpret_cast<fftw_complex*> (complexResponse.memptr());
    fftw_plan plan2 = fftw_plan_dft_2d(complexResponse.n_cols,complexResponse.n_rows,
                                      in,  out, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(plan2);

    complexResponse = Functions::fftShift(complexResponse);
    m_response = real(complexResponse);

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




//void Neuron::computeResponse(double t)
//{

//    complexResponse.set_real(m_stim->spatial());
//    complexResponse.set_imag(m_stim->spatial()*0);
//    complexResponse = Functions::fftShift(complexResponse);


//    fftw_complex* in = reinterpret_cast<fftw_complex*> (complexResponse.memptr());
//    fftw_complex* out = reinterpret_cast<fftw_complex*> (complexResponse.memptr());

//    fftw_plan plan = fftw_plan_dft_2d(complexResponse.n_cols,complexResponse.n_rows,
//                                       in,  out, FFTW_FORWARD, FFTW_ESTIMATE);

//    fftw_execute(plan);
//    complexResponse /=m_freqMesh.n_elem*m_freqMesh.n_elem;

//    fftw_plan plan1 = fftw_plan_dft_2d(complexResponse.n_cols,complexResponse.n_rows,
//                                       in,  out, FFTW_BACKWARD, FFTW_ESTIMATE);

//    fftw_execute(plan1);

//    m_response =real((complexResponse));
//    m_response = Functions::fftShift(m_response);

//}





//void Neuron::computeResponse(double t)
//{

//    mat stim = 0*m_response;
//    m_response = 0*m_response;
//    double *w = new double [int(m_domain[2])];
//    double *x = new double [int(m_domain[2])];
//    gauleg(m_domain[0], m_domain[1], x, w, m_domain[2]);

//    for(int i = 0; i < int(m_spatialMesh.n_elem); i++){
//        for(int j = 0; j < int(m_spatialMesh.n_elem); j++){

//            stim(i,j) = m_stim->real({m_spatialMesh[i], m_spatialMesh[j]}, t);
//            for(int m = 0; m < int(m_domain[2]); m++){
//                for(int n = 0; n < int(m_domain[2]); n++){

//                    double Scomplex = m_stim->complex({x[m], x[n]}, m_stim->w());
//                    double Gcomplex = impulseResponseComplex({x[m], x[n]}, m_stim->w());

//                    double dGr = 1./8./(PI*PI*PI) * Gcomplex * w[m] * w[n] *
//                            cos(m_spatialMesh[i]*x[m]+ m_spatialMesh[j]*x[n] - m_stim->w() * t);
//                    m_response(i,j) += dGr * Scomplex;

//                }
//            }
//        }
//    }
//    m_stim->setReal(stim);

//}

void Neuron::computeResponseComplex(double w)
{
    cout << "computeResponseComplex: Not implemented!" << endl;
    for(int i = 0; i < int(m_spatialMesh.n_elem); i++){
        for(int j = 0; j < int(m_spatialMesh.n_elem); j++){

            double Gcomplex = impulseResponseComplex({m_spatialMesh[i], m_spatialMesh[j]}, w);
            double Scomplex = m_stim->frequency({m_spatialMesh[i], m_spatialMesh[j]}, w);

            m_impulseResponseComplex(i,j) = Gcomplex;
            m_responseComplex(i,j) = Gcomplex*Scomplex;
        }
    }

}

void Neuron::computeImpulseResponse(double t)
{
    m_impulseResponse  = 0* m_impulseResponse;
    double *w = new double [int(m_domain[2])];
    double *x = new double [int(m_domain[2])];
    gauleg(m_domain[0], m_domain[1], x, w, m_domain[2]);

    for(int i = 0; i < int(m_spatialMesh.n_elem); i++){
        for(int j = 0; j < int(m_spatialMesh.n_elem); j++){

            for(int m = 0; m < int(m_domain[2]); m++){
                for(int n = 0; n < int(m_domain[2]); n++){
                    for(int o = 0; o < int(m_domain[2]); o++){
                        m_impulseResponse(i, j) += 1./8./(PI*PI*PI) *
                                impulseResponseComplex({x[m], x[n]},x[o])
                                * w[m] * w[n] * w[o]
                                * cos(m_spatialMesh[i]*x[m]+ m_spatialMesh[j]*x[n] - x[o]* t);
                    }

                }
            }
        }
    }
}

void Neuron::computeImpulseResponseComplex(double w)
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

mat Neuron::responseComplex() const
{
    return m_responseComplex;
}

mat Neuron::impulseResponseComplex() const
{
    return m_impulseResponseComplex;
}










