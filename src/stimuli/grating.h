#ifndef GRATING_H
#define GRATING_H

#include "stimuli.h"


class Grating : public Stimuli
{
public:
    Grating(const Config *cfg, Integrator integrator);
    ~Grating();

    // Stimuli interface
private:
    double m_contrast = 0.0;

    double valueAtPoint(vec2 rVec, double t);
    double fourierTransformAtFrequency(vec2 k, double w);

};

#endif // GRATING_H
