#include "neuron.h"

Neuron::Neuron(const Config *cfg, Stimuli *stim)
    : m_stim(stim)
{
    const Setting & root = cfg->getRoot();
    m_nPoints = root["spatialDomainSettings"]["nPoints"];
    m_nSteps = root["dynamicSettings"]["nSteps"];

    m_response = zeros(m_nPoints, m_nPoints, m_nSteps);
    m_impulseResponse = zeros(m_nPoints, m_nPoints, m_nSteps);
    m_responseFT = zeros<cx_cube>(m_nPoints, m_nPoints, m_nSteps);
    m_impulseResponseFT = zeros<cx_cube>(m_nPoints, m_nPoints, m_nSteps);

    //Spatial Mesh
    m_spatialMesh = linspace(-0.5, 0.5, m_nPoints);
    double dr = m_spatialMesh(1) - m_spatialMesh(0);
    double N_2 = ceil(m_nPoints/2.);
    double df = 1./dr/m_nPoints;
    double fs = 1./dr;

    m_spatialFreqs = linspace(-N_2*df, (m_nPoints - 1. - N_2)*df, m_nPoints);
    m_spatialFreqs*= 2*PI;


    //Temporal Mesh
    m_temporalMesh = linspace(-0.5, 0.5, m_nSteps);
    double dt = m_temporalMesh(1) - m_temporalMesh(0);
    double Nt_2 = ceil(m_nSteps/2.);
    double df_t = 1./dt/m_nSteps;
    double fs_t = 1./dt;

    m_temporalFreqs = linspace(-Nt_2*df_t, (m_nSteps - 1. - Nt_2)*df_t, m_nSteps);
    m_temporalFreqs*= 2*PI;
}

Neuron::~Neuron()
{

}


void Neuron::computeResponse()
{
    computeImpulseResponseFT();

    m_responseFT = /*m_impulseResponseFT % */m_stim->fourierTransform();
    m_responseFT = Functions::fftShift3d(m_responseFT);

    int size[3] = {m_nPoints, m_nPoints , m_nSteps};
    fftw_complex* in = reinterpret_cast<fftw_complex*> (m_responseFT.memptr());
    fftw_complex* out = reinterpret_cast<fftw_complex*> (m_responseFT.memptr());
    fftw_plan plan = fftw_plan_dft(3, size, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(plan);
    fftw_destroy_plan(plan);

    m_responseFT = Functions::fftShift3d(m_responseFT);
    m_response = real(m_responseFT);

    for(int k = 0; k < int(m_response.n_slices); k++){
            m_response.slice(k) = flipud(m_response.slice(k));
    }
}


void Neuron::computeImpulseResponse()
{
    computeImpulseResponseFT();

    m_impulseResponseFT = Functions::fftShift3d(m_impulseResponseFT);
    int size[3] = {m_nPoints, m_nPoints , m_nSteps};
    fftw_complex* in = reinterpret_cast<fftw_complex*> (m_impulseResponseFT.memptr());
    fftw_complex* out = reinterpret_cast<fftw_complex*> (m_impulseResponseFT.memptr());
    fftw_plan plan = fftw_plan_dft(3, size, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(plan);
    fftw_destroy_plan(plan);

    m_impulseResponseFT = Functions::fftShift3d(m_impulseResponseFT);
    m_impulseResponse = real(m_impulseResponseFT);
}


void Neuron::computeImpulseResponseFT()
{
    for(int k = 0; k < m_nSteps; k++){
        for(int i = 0; i < m_nPoints; i++){
            for(int j = 0; j < m_nPoints; j++){
                m_impulseResponseFT(i,j,k) =
                        impulseResponseFT({m_spatialFreqs[i],m_spatialFreqs[j]},
                                          m_temporalFreqs[k]);
            }
        }
    }
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

cube Neuron::response() const
{
    return m_response;
}

cube Neuron::impulseResponse() const
{
    return m_impulseResponse;
}

cx_cube Neuron::responseFT() const
{
    return m_responseFT;
}

cx_cube Neuron::impulseResponseFT() const
{
    return m_impulseResponseFT;
}










