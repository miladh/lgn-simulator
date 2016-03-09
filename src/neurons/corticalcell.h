#ifndef CORTICALCELL_H
#define CORTICALCELL_H

#include "neuron.h"

namespace lgnSimulator {
class CorticalCell : public Neuron
{

public:
    CorticalCell(const Integrator &integrator, double backgroundResponse = 0);
    ~CorticalCell();

    void addRelayCell(Neuron* const neuron, const Kernel &kernell);
    const Kernel *relayInputKernel() const;

    // Neuron interface
    virtual void computeImpulseResponseFourierTransform();

private:
    Input* m_relayInput = nullptr;
    void computeNeededcubes() const;

};
}
#endif // CORTICALCELL_H
