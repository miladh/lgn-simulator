#ifndef GRATING_H
#define GRATING_H

#include "stimuli.h"


class Grating : public Stimulus
{
public:
    Grating(Integrator integrator, vec2 kd, double wd, double contrast);
    ~Grating();

    // Stimuli interface
private:
    double m_contrast = 0.0;
    vec2 m_k = {0,0};
    double m_w = 0;

    double valueAtPoint(vec2 rVec, double t);
    double fourierTransformAtFrequency(vec2 k, double w);

};

Grating *createGratingStimulus(Integrator integrator, const Config *cfg);

#endif // GRATING_H
