#ifndef GRATING_H
#define GRATING_H

#include "stimuli.h"


class Grating : public Stimulus
{
public:
    Grating(Integrator *integrator, vec2 kd, double wd, double contrast);
    ~Grating();

    // Stimulus interface
public:
    virtual void computeSpatiotemporal();
    virtual void computeFourierTransform();
private:
    double valueAtPoint(vec2 rVec, double t);
    double fourierTransformAtFrequency(vec2 kVec, double w);

private:
    vec2 m_k = {0,0};
    double m_w = 0;
    double m_contrast = 0.0;

};

Grating createGratingStimulus(Integrator *integrator, const Config *cfg);

#endif // GRATING_H
