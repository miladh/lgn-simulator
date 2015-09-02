#include "relaycell.h"

RelayCell::RelayCell(const Config *cfg, Stimuli *stim)
    : Neuron(cfg, stim)

{
    m_cellType = "relay";
}

RelayCell::~RelayCell()
{

}


double RelayCell::impulseResponseFT(vec2 kVec, double w)
{

    double G = 0;
    double I = 0;
    double C = 0;

    for (const Input g : m_ganglionCells){
        Neuron *ganglionCell = g.neuron;
        G += g.spatialKernel->complex(kVec)
                * g.temporalKernel->complex(w)
                * ganglionCell->impulseResponseFT(kVec,w);
    }



    for (const Input i : m_interNeurons){
        Neuron *interneuron = i.neuron;
        double Kri = i.spatialKernel->complex(kVec)
                * i.temporalKernel->complex(w);


        for (const Input g : interneuron->ganglionCells()){
            Neuron *ganglionCell = g.neuron;
            I += g.spatialKernel->complex(kVec)
                    * g.temporalKernel->complex(w)
                    * ganglionCell->impulseResponseFT(kVec,w);
        }
        I*= Kri;
    }


    for (const Input c : m_corticalNeurons){
        Neuron *corticalCell = c.neuron;
        double Krc = c.spatialKernel->complex(kVec)
                * c.temporalKernel->complex(w);

        for (const Input r : corticalCell->relayCells()){
            C += r.spatialKernel->complex(kVec)
                    * r.temporalKernel->complex(w);
        }
        C*= Krc;
    }

    double Gr = (G + I)/(1 - C);
    return Gr;
}

