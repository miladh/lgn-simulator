#ifndef DOGSTIM_H
#define DOGSTIM_H

#include "stimuli.h"
#include "../spatialKernels/dog.h"


class DOGstim : public Stimuli
{
public:
    DOGstim(const Config *cfg, Integrator integrator);
    ~DOGstim();

    // Stimuli interface
public:

    DOG m_dog = DOG(1.0, 0.25, 0.85, 0.83);

//    DOG m_dog = DOG(1.0, 0.25, 0.0, 0.83);

private:
    double valueAtPoint(vec2 rVec, double t);
    double fourierTransformAtFrequency(vec2 k, double w);

};

#endif // DOGSTIM_H
