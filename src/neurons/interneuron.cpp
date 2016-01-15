#include "interneuron.h"


using namespace edog;


Interneuron::Interneuron(Integrator *integrator)
    : Neuron(integrator)
{
    m_cellType = "interneuron";
}

Interneuron::~Interneuron()
{

}

void Interneuron::computeImpulseResponseFourierTransform()
{
    computeNeededcubes();
    impulseResponseFourierTransformComputed = true;

    for(int k = 0; k < int(m_impulseResponseFT.n_slices); k++){
        double w = -m_temporalFreqs[k];

        for(int i = 0; i < int(m_impulseResponseFT.n_rows); i++){
            for(int j = 0; j < int(m_impulseResponseFT.n_cols); j++){
                vec2 kVec= {m_spatialFreqs[i], m_spatialFreqs[j]};

                for (const Input g : m_ganglionCells){
                    Neuron *ganglionCell = g.neuron;
                    m_impulseResponseFT(i,j,k)
                            += g.spatialKernel->fourierTransform(kVec)
                            * g.temporalKernel->fourierTransform(w)
                            * ganglionCell->impulseResponseFourierTransform()(i,j,k);
                }

                for (const Input c : m_corticalNeurons){
                    Neuron *corticalCell = c.neuron;
                    m_impulseResponseFT(i,j,k)
                            += c.spatialKernel->fourierTransform(kVec)
                            * c.temporalKernel->fourierTransform(w)
                            * corticalCell->impulseResponseFourierTransform()(i,j,k);
                }
            }
        }
    }
}

void Interneuron::computeNeededcubes()
{
    for (const Input g : m_ganglionCells){
        Neuron *ganglionCell = g.neuron;
        if(!ganglionCell->isImpulseResponseFourierTransformComputed()){
            ganglionCell->computeImpulseResponseFourierTransform();
        }
    }

    for (const Input c : m_corticalNeurons){
        Neuron *corticalCell = c.neuron;
        if(!corticalCell->isImpulseResponseFourierTransformComputed()){
            corticalCell->computeImpulseResponseFourierTransform();
        }
    }
}

