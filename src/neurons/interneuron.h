#ifndef INTERNEURONS_H
#define INTERNEURONS_H

#include "neuron.h"


class Interneuron : public Neuron
{
    friend class RelayCell;

public:
    Interneuron(Integrator *integrator);
    ~Interneuron();

    // Neuron interface
    virtual void computeImpulseResponseFourierTransform();
};

#endif // INTERNEURONS_H
