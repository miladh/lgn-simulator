#ifndef RELAYCELL_H
#define RELAYCELL_H

#include "neuron.h"

class RelayCell : public Neuron
{
public:
    RelayCell(const Config * cfg);
    ~RelayCell();

    // Neuron interface
    void computeResponse(double t);

private:
    double transferFunctionComplex(vec2 kVec, double w);
};

#endif // RELAYCELL_H
