#include "neuron.h"

Neuron::Neuron(Integrator *integrator, StaticNonlinearity *staticNonlinearity)
    : m_integrator(integrator)
    , m_staticNonlinearity(staticNonlinearity)
{
    int nPointsTemporal = integrator->nPointsTemporal();
    int nPointsSpatial = integrator->nPointsSpatial();

    m_impulseResponseFT=zeros<cx_cube>(nPointsSpatial,nPointsSpatial, nPointsTemporal);

    //Temporal Mesh
    m_timeVec = integrator->timeVec();
    m_temporalFreqs = integrator->temporalFreqVec();

    //Spatial Mesh
    m_coordinateVec = integrator->coordinateVec();
    m_spatialFreqs =integrator->spatialFreqVec();
}

Neuron::~Neuron()
{

}


void Neuron::computeResponse(Stimulus *stimulus)
{
    if(!impulseResponseFourierTransformComputed){
        computeImpulseResponseFourierTransform();
    }
    m_responseFT = m_impulseResponseFT % stimulus->fourierTransform();
    m_response = real(m_integrator->backwardFFT(m_responseFT));

    if(m_staticNonlinearity != nullptr){
        m_staticNonlinearity->applyStaticNonlinearity(&m_response);
    }


}

void Neuron::computeImpulseResponse()
{
    if(!impulseResponseFourierTransformComputed){
        computeImpulseResponseFourierTransform();
    }
    m_impulseResponse = real(m_integrator->backwardFFT(m_impulseResponseFT));

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



bool Neuron::isImpulseResponseFourierTransformComputed() const
{
    return impulseResponseFourierTransformComputed;
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


const cube& Neuron::response() const
{
    return m_response;
}

const cube& Neuron::impulseResponse() const
{
    return m_impulseResponse;
}

const cx_cube& Neuron::responseFT() const
{
    return m_responseFT;
}

const cx_cube& Neuron::impulseResponseFourierTransform() const
{
    return m_impulseResponseFT;
}


void Neuron::clearResponse()
{
    m_response.clear();
    m_responseFT.clear();
}


void Neuron::clearImpulseResponse()
{
    m_impulseResponse.clear();

}







