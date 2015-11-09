#include "corticalcell.h"

CorticalCell::CorticalCell(Integrator* integrator)
    : Neuron(integrator)
{
        m_cellType = "cortical";

}

CorticalCell::~CorticalCell()
{

}

complex<double> CorticalCell::impulseResponseFourierTransformAtFrequency(vec2 kVec, double w)
{
    complex<double> R = 0;
    for (const Input r : m_relayCells){
        Neuron *relayCell = r.neuron;
        R += r.spatialKernel->fourierTransform(kVec)
                * r.temporalKernel->fourierTransform(w)
                * relayCell->impulseResponseFourierTransformAtFrequency(kVec,w);
    }

    return R;
}


