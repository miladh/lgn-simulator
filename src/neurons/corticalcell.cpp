#include "corticalcell.h"

using namespace lgnSimulator;


CorticalCell::CorticalCell(Integrator* const integrator, double backgroundResponse,
                           StaticNonlinearity * const staticNonlinearity)
    : Neuron(integrator, backgroundResponse, "cortical", staticNonlinearity)
{
}

CorticalCell::~CorticalCell()
{
}

void CorticalCell::computeImpulseResponseFourierTransform()
{
    computeNeededcubes();
    impulseResponseFourierTransformComputed = true;

    for(int k = 0; k < int(m_impulseResponseFT.n_slices); k++){
        double w = m_temporalFreqs[k];

        for(int i = 0; i < int(m_impulseResponseFT.n_rows); i++){
            for(int j = 0; j < int(m_impulseResponseFT.n_cols); j++){
                vec2 kVec= {m_spatialFreqs[i], m_spatialFreqs[j]};

                m_impulseResponseFT(i,j,k)
                        = m_relayInput->kernel.fourierTransform(kVec,w)
                        * m_relayInput->
                        neuron->impulseResponseFourierTransform()(i,j,k);
            }
        }
    }
}

void CorticalCell::computeNeededcubes() const
{
    if(!m_relayInput->neuron->isImpulseResponseFourierTransformComputed()){
        m_relayInput->neuron->computeImpulseResponseFourierTransform();
    }
}



void CorticalCell::addRelayCell(Neuron* const neuron, const Kernel &kernel)
{
    if (neuron->type() == "relay") {
        m_relayInput = new Input{neuron, kernel};
    }else{
        throw overflow_error("wrong cell type in addRelayCell(): "
                             + neuron->type());
    }
}

const Kernel* CorticalCell::relayInputKernel() const
{
    return &m_relayInput->kernel;
}
