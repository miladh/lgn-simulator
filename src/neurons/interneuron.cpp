#include "interneuron.h"



Interneuron::Interneuron(const Config *cfg, Stimuli *stim, Integrator integrator)
    : Neuron(cfg, stim, integrator)
{
    m_cellType = "interneuron";
}

Interneuron::~Interneuron()
{

}



double Interneuron::impulseResponseFT(vec2 kVec, double w)
{


    double G = 0;
    double C = 0;

    for (const Input g : m_ganglionCells){
        Neuron *ganglionCell = g.neuron;
        G += g.spatialKernel->fourierTransform(kVec)
                * g.temporalKernel->complex(w)
                * ganglionCell->impulseResponseFT(kVec,w);
    }

    for (const Input c : m_corticalNeurons){
        Neuron *corticalCell = c.neuron;
        double Krc = c.spatialKernel->fourierTransform(kVec)
                * c.temporalKernel->complex(w);

        for (const Input r : corticalCell->relayCells()){
            C += r.spatialKernel->fourierTransform(kVec)
                    * r.temporalKernel->complex(w);
        }
        C*= Krc;
    }

    double Gr = (G + C);
    return Gr;
}
