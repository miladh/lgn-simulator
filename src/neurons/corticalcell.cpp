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
                        = m_relayInput.at(0).kernel.fourierTransform(kVec,w)
                        * m_relayInput.at(0).
                        neuron->impulseResponseFourierTransform()(i,j,k);
            }
        }
    }
}

void CorticalCell::computeNeededcubes() const
{
    if(!m_relayInput.at(0).neuron->isImpulseResponseFourierTransformComputed()){
        m_relayInput.at(0).neuron->computeImpulseResponseFourierTransform();
    }
}



void CorticalCell::addRelayCell(Neuron* const neuron, const Kernel &kernel)
{
    if(m_relayInput.size()>1){
        throw length_error("cortical cell already has relay cell input");

    }else if (neuron->type() == "relay") {
        m_relayInput.emplace_back(Input{neuron, kernel});

    }else{
        throw overflow_error("wrong cell type in addRelayCell(): "
                             + neuron->type());
    }
}

const Kernel* CorticalCell::relayInputKernel() const
{
    return &m_relayInput.at(0).kernel;
}
