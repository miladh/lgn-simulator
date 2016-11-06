#include "neuron.h"

/*!
 \class lgnSimulator::Neuron
 \inmodule lgnSimulator
 \brief Virtual class for neurons.
 */

using namespace lgnSimulator;



Neuron::Neuron(Integrator* const integrator,
               const double backgroundResponse,
               const string type,
               StaticNonlinearity * const staticNonlinearity)
    : m_integrator(integrator)
    , m_staticNonlinearity(staticNonlinearity)
    , m_backgroundResponse(backgroundResponse)
    , m_type(type)
{
    int nPointsTemporal = integrator->nPointsTemporal();
    int nPointsSpatial = integrator->nPointsSpatial();

    m_impulseResponseFT=zeros<cx_cube>(nPointsSpatial, nPointsSpatial, nPointsTemporal);

    //Temporal Mesh
    m_timeVec = integrator->timeVec();
    m_temporalFreqs = integrator->temporalFreqVec();

    //Spatial Mesh
    m_spatialVec = integrator->spatialVec();
    m_spatialFreqs =integrator->spatialFreqVec();
}

Neuron::~Neuron()
{

}

void Neuron::computeResponse(const Stimulus *const stimulus,
                             bool recomputeFourierTransfrom)
{
    if(!responseFourierTransformComputed || recomputeFourierTransfrom){
        computeResponseFourierTransform(stimulus);
    }

    m_response = m_integrator->backwardFFT(m_responseFT);

    if(m_staticNonlinearity != nullptr){
        m_staticNonlinearity->applyStaticNonlinearity(&m_response);
    }
}

void Neuron::computeResponseFourierTransform(const Stimulus * const stimulus)
{
    responseFourierTransformComputed=true;

    if(!impulseResponseFourierTransformComputed){
        computeImpulseResponseFourierTransform();
    }

    m_responseFT = m_impulseResponseFT % stimulus->fourierTransform();

    if(m_backgroundResponse!=0.0){ //add DC contribution from bck activity
        m_responseFT(0,0,0) += 8*core::pi*core::pi*core::pi * m_backgroundResponse
                /m_integrator->spatialFreqResolution()
                /m_integrator->spatialFreqResolution()
                /m_integrator->temporalFreqResolution();
    }
}

void Neuron::computeImpulseResponse()
{
    if(!impulseResponseFourierTransformComputed){
        computeImpulseResponseFourierTransform();
    }
    m_impulseResponse = m_integrator->backwardFFT(m_impulseResponseFT);

}


bool Neuron::isImpulseResponseFourierTransformComputed() const
{
    return impulseResponseFourierTransformComputed;
}


const string Neuron::type() const
{
    return m_type;
}


const cube& Neuron::response() const
{
    return m_response;
}

const cube& Neuron::impulseResponse() const
{
    return m_impulseResponse;
}

const cx_cube& Neuron::responseFourierTransform() const
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
}

void Neuron::clearResponseFourierTransform()
{
    responseFourierTransformComputed=false;
    m_responseFT.clear();
}


void Neuron::clearImpulseResponse()
{
    m_impulseResponse.clear();

}







