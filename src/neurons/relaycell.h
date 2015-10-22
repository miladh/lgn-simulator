#ifndef RELAYCELL_H
#define RELAYCELL_H

#include "neuron.h"

class RelayCell : public Neuron
{
public:
    RelayCell(Integrator *integrator);
    ~RelayCell();

private:
    // Neuron interface
    double impulseResponseFourierTransformAtFrequency(vec2 kVec, double w);
};

#endif // RELAYCELL_H
