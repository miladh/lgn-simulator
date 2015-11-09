#include "interneuron.h"



Interneuron::Interneuron(Integrator *integrator)
    : Neuron(integrator)
{
    m_cellType = "interneuron";
}

Interneuron::~Interneuron()
{

}



complex<double> Interneuron::impulseResponseFourierTransformAtFrequency(vec2 kVec, double w)
{


    complex<double> G = 0;
    complex<double> C = 0;

    for (const Input g : m_ganglionCells){
        Neuron *ganglionCell = g.neuron;
        G += g.spatialKernel->fourierTransform(kVec)
                * g.temporalKernel->fourierTransform(w)
                * ganglionCell->impulseResponseFourierTransformAtFrequency(kVec,w);
    }

    for (const Input c : m_corticalNeurons){
        Neuron *corticalCell = c.neuron;
        complex<double> Kic = c.spatialKernel->fourierTransform(kVec)
                * c.temporalKernel->fourierTransform(w);

        for (const Input r : corticalCell->relayCells()){
            Neuron *relayCell = r.neuron;
            C += r.spatialKernel->fourierTransform(kVec)
                    * r.temporalKernel->fourierTransform(w)
                    * relayCell->impulseResponseFourierTransformAtFrequency(kVec,w);
        }
        C*= Kic;
    }

    complex<double> Gr = (G + C);
    return Gr;
}
