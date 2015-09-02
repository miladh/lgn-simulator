#include "corticalcell.h"

CorticalCell::CorticalCell(const Config * cfg, Stimuli *stim)
    : Neuron(cfg, stim)
{
        m_cellType = "cortical";

}

CorticalCell::~CorticalCell()
{

}

double CorticalCell::impulseResponseFT(vec2 kVec, double w)
{
    double R = 0;
    for (const Input r : m_relayCells){
        Neuron *relayCell = r.neuron;
        R += r.spatialKernel->complex(kVec)
                * r.temporalKernel->complex(w)
                * relayCell->impulseResponseFT(kVec,w);
    }

    return R;
}


