#ifndef CORTICALCELL_H
#define CORTICALCELL_H

#include "neuron.h"

namespace lgnSimulator {
class CorticalCell : public Neuron
{
    friend class RelayCell;

public:
    CorticalCell(const Integrator &integrator, double backgroundResponse = 0);
    ~CorticalCell();

    // Neuron interface
    virtual void computeImpulseResponseFourierTransform();

private:
    void computeNeededcubes();
};
}
#endif // CORTICALCELL_H
