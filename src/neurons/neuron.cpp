#include "neuron.h"

Neuron::Neuron(const Config *cfg, Stimuli *stim, Integrator integrator)
    : m_stim(stim)
    , m_integrator(integrator)
{
    const Setting & root = cfg->getRoot();
    int nPointsTemporal = integrator.nPointsTemporal();
    int nPointsSpatial = integrator.nPointsSpatial();

    m_response = zeros(nPointsSpatial, nPointsSpatial, nPointsTemporal);
    m_impulseResponse = zeros(nPointsSpatial, nPointsSpatial, nPointsTemporal);
    m_responseFT = zeros<cx_cube>(nPointsSpatial, nPointsSpatial, nPointsTemporal);
    m_impulseResponseFT=zeros<cx_cube>(nPointsSpatial,nPointsSpatial, nPointsTemporal);

    //Temporal Mesh
    timeVec = integrator.timeVec();
    m_temporalFreqs = integrator.temporalFreqVec();

    //Spatial Mesh
    m_coordinateVec = integrator.coordinateVec();
    m_spatialFreqs =integrator.spatialFreqVec();

}

Neuron::~Neuron()
{

}


void Neuron::computeResponse()
{
    computeImpulseResponseFT();

    m_responseFT = m_impulseResponseFT % m_stim->fourierTransform();
    m_responseFT = m_integrator.integrate(m_responseFT);
    m_responseFT = FFTHelper::fftShift(m_responseFT);

    m_response = real(m_responseFT);
}


void Neuron::computeImpulseResponse()
{
    computeImpulseResponseFT();

    m_impulseResponseFT = m_integrator.integrate(m_impulseResponseFT);
    m_impulseResponseFT = FFTHelper::fftShift(m_impulseResponseFT);

    m_impulseResponse = real(m_impulseResponseFT);
}


void Neuron::computeImpulseResponseFT()
{
    for(int k = 0; k < int(m_impulseResponseFT.n_slices); k++){
        for(int i = 0; i < int(m_impulseResponseFT.n_rows); i++){
            for(int j = 0; j < int(m_impulseResponseFT.n_cols); j++){
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










