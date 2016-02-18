#ifndef INTERNEURONS_H
#define INTERNEURONS_H

#include "neuron.h"

namespace lgnSimulator {
class Interneuron : public Neuron
{
    friend class RelayCell;

public:
    Interneuron(const Integrator &integrator, double backgroundResponse= 0);
    ~Interneuron();

    // Neuron interface
    virtual void computeImpulseResponseFourierTransform();

private:
    void computeNeededcubes();
};
}
#endif // INTERNEURONS_H
