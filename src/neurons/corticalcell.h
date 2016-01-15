#ifndef CORTICALCELL_H
#define CORTICALCELL_H

#include "neuron.h"

namespace edog {
class CorticalCell : public Neuron
{
    friend class RelayCell;

public:
    CorticalCell(Integrator *integrator);
    ~CorticalCell();

    // Neuron interface
    virtual void computeImpulseResponseFourierTransform();

private:
    void computeNeededcubes();
};
}
#endif // CORTICALCELL_H
