#ifndef INTERNEURONS_H
#define INTERNEURONS_H

#include "neuron.h"



class Interneuron : public Neuron
{
public:
    Interneuron(Integrator *integrator);
    ~Interneuron();

    // Neuron interface
public:
    virtual double impulseResponseFT(vec2 kVec, double w);
};

#endif // INTERNEURONS_H
