#ifndef INTERNEURONS_H
#define INTERNEURONS_H

#include "neuron.h"

namespace lgnSimulator {
class Interneuron : public Neuron
{
    friend class RelayCell;

public:
    Interneuron(Integrator *integrator);
    ~Interneuron();

    // Neuron interface
    virtual void computeImpulseResponseFourierTransform();

private:
    void computeNeededcubes();
};
}
#endif // INTERNEURONS_H
