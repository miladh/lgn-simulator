#ifndef CORTICALCELL_H
#define CORTICALCELL_H

#include "neuron.h"

class CorticalCell : public Neuron
{
public:
    CorticalCell(Integrator *integrator);
    ~CorticalCell();

protected:

    // Neuron interface
    double impulseResponseFT(vec2 kVec, double w);
};

#endif // CORTICALCELL_H
