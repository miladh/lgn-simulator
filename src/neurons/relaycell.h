#ifndef RELAYCELL_H
#define RELAYCELL_H

#include "neuron.h"
namespace lgnSimulator {
class RelayCell : public Neuron
{
public:
    RelayCell(Integrator *integrator, double backgroundResponse= 0);
    ~RelayCell();

    // Neuron interface
    virtual void computeImpulseResponseFourierTransform();

private:
    complex<double> impulseResponseFourierTransformAtFrequency(int idx, int jdx, int kdx);
    void computeNeededcubes();
};
}
#endif // RELAYCELL_H
