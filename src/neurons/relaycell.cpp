#include "relaycell.h"

RelayCell::RelayCell(const Config *cfg, Stimuli *stim, Integrator integrator)
    : Neuron(cfg, stim, integrator)

{
    m_cellType = "relay";
}

RelayCell::~RelayCell()
{

}


double RelayCell::impulseResponseFT(vec2 kVec, double w)
{

    double G = 0;
    double Iff = 0;
    double Ifb = 0;
    double C = 0;

    for (const Input g : m_ganglionCells){
        Neuron *ganglionCell = g.neuron;
        G += g.spatialKernel->fourierTransform(kVec)
                * g.temporalKernel->fourierTransform(w)
                * ganglionCell->impulseResponseFT(kVec,w);
    }



    for (const Input i : m_interNeurons){
        Neuron *interneuron = i.neuron;
        double Kri = i.spatialKernel->fourierTransform(kVec)
                * i.temporalKernel->fourierTransform(w);


        //Feedforward term
        for (const Input g : interneuron->ganglionCells()){
            Neuron *ganglionCell = g.neuron;
            Iff += g.spatialKernel->fourierTransform(kVec)
                    * g.temporalKernel->fourierTransform(w)
                    * ganglionCell->impulseResponseFT(kVec,w);
        }
        Iff*= Kri;

        //Feedback term
        for (const Input c : interneuron->corticalNeurons()){
            Neuron *corticalCell = c.neuron;
            double Kic = c.spatialKernel->fourierTransform(kVec)
                    * c.temporalKernel->fourierTransform(w);

            // NOTE: ONLY ONE RELAY CELL!!!
            for (const Input r : corticalCell->relayCells()){
                double Kcr = r.spatialKernel->fourierTransform(kVec)
                        * r.temporalKernel->fourierTransform(w);
                Ifb += Kri*Kic*Kcr;
            }

        }

    }



    for (const Input c : m_corticalNeurons){
        Neuron *corticalCell = c.neuron;
        double Krc = c.spatialKernel->fourierTransform(kVec)
                * c.temporalKernel->fourierTransform(w);

        for (const Input r : corticalCell->relayCells()){
            C += r.spatialKernel->fourierTransform(kVec)
                    * r.temporalKernel->fourierTransform(w);
        }
        C*= Krc;
    }


    double Gr = (G + Iff)/(1 - Ifb - C);
    return Gr;
}

