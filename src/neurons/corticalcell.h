#ifndef CORTICALCELL_H
#define CORTICALCELL_H

#include "neuron.h"

class CorticalCell : public Neuron
{
public:
    CorticalCell(const Config *cfg, Stimuli *stim);
    ~CorticalCell();

protected:

    // Neuron interface
    double impulseResponseComplex(vec2 kVec, double w);
};

#endif // CORTICALCELL_H
