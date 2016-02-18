#include "neuron.h"

using namespace lgnSimulator;

/*!
 * \class Neuron
 * \inmodule lgnSimulator
 * \brief Virtual class for neurons.
 *
 */


Neuron::Neuron(const Integrator& integrator,
               const double backgroundResponse,
               StaticNonlinearity * const staticNonlinearity)
    : m_integrator(integrator)
    , m_staticNonlinearity(staticNonlinearity)
    , m_backgroundResponse(backgroundResponse)
{
    int nPointsTemporal = integrator.nPointsTemporal();
    int nPointsSpatial = integrator.nPointsSpatial();

    m_impulseResponseFT=zeros<cx_cube>(nPointsSpatial,nPointsSpatial, nPointsTemporal);

    //Temporal Mesh
    m_timeVec = integrator.timeVec();
    m_temporalFreqs = integrator.temporalFreqVec();

    //Spatial Mesh
    m_spatialVec = integrator.spatialVec();
    m_spatialFreqs =integrator.spatialFreqVec();
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
    if(m_backgroundResponse!=0){ //add DC contribution from bck activity
        m_responseFT(0,0,0) += 8*core::pi*core::pi*core::pi * m_backgroundResponse
            /m_integrator.spatialFreqResolution()
            /m_integrator.spatialFreqResolution()
            /m_integrator.temporalFreqResolution();
    }

    m_response = real(m_integrator.backwardFFT(m_responseFT));

    if(m_staticNonlinearity != nullptr){
        m_staticNonlinearity->applyStaticNonlinearity(&m_response);
    }


}

void Neuron::computeImpulseResponse()
{
    if(!impulseResponseFourierTransformComputed){
        computeImpulseResponseFourierTransform();
    }
    m_impulseResponse = real(m_integrator.backwardFFT(m_impulseResponseFT));

}

void Neuron::addGanglionCell(Neuron* const neuron, const Kernel &kernel)
{
    m_ganglionCells.emplace_back(Input{neuron, kernel});
}

void Neuron::addRelayCell(Neuron* const neuron, const Kernel &kernel)
{

    m_relayCells.emplace_back(Input{neuron, kernel});
}

void Neuron::addInterNeuron(Neuron* const neuron, const Kernel &kernel)
{
    m_interNeurons.emplace_back(Input{neuron, kernel});
}

void Neuron::addCorticalNeuron(Neuron* const neuron, const Kernel &kernel)
{
    m_corticalNeurons.emplace_back(Input{neuron, kernel});
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







