#include "corticalcell.h"

using namespace lgnSimulator;


CorticalCell::CorticalCell(const Integrator& integrator, double backgroundResponse)
    : Neuron(integrator, backgroundResponse)
{
    m_cellType = "cortical";

}

CorticalCell::~CorticalCell()
{
}

void CorticalCell::computeImpulseResponseFourierTransform()
{
    computeNeededcubes();
    impulseResponseFourierTransformComputed = true;

    for(int k = 0; k < int(m_impulseResponseFT.n_slices); k++){
        double w = -m_temporalFreqs[k];

        for(int i = 0; i < int(m_impulseResponseFT.n_rows); i++){
            for(int j = 0; j < int(m_impulseResponseFT.n_cols); j++){
                vec2 kVec= {m_spatialFreqs[i], m_spatialFreqs[j]};

                for (const Input r : m_relayCells){
                    Neuron *relayCell = r.neuron;
                    m_impulseResponseFT(i,j,k)
                            += r.kernel.fourierTransform(kVec,w)
                            * relayCell->impulseResponseFourierTransform()(i,j,k);
                }
            }
        }
    }
}

void CorticalCell::computeNeededcubes()
{
    for (const Input r : m_relayCells){
        Neuron *relayCell = r.neuron;
        if(!relayCell->isImpulseResponseFourierTransformComputed()){
            relayCell->computeImpulseResponseFourierTransform();
        }
    }
}


