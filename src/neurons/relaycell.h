#ifndef RELAYCELL_H
#define RELAYCELL_H

#include "neuron.h"

class RelayCell : public Neuron
{
public:
    RelayCell(const Config * cfg, Stimuli *stim);
    ~RelayCell();

    // Neuron interface
    void computeResponse(double t);
    double impulseResponseComplex(vec2 kVec, double w);
    void computeResponseComplex(double w);
};

#endif // RELAYCELL_H
