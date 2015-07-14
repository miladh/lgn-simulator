#ifndef CORTICALCELL_H
#define CORTICALCELL_H

#include "neuron.h"

class CorticalCell : public Neuron
{
public:
    CorticalCell(const Config *cfg);
    ~CorticalCell();

    // Neuron interface
    void computeResponse();
};

#endif // CORTICALCELL_H
