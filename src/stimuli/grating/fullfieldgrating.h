#ifndef FULLFIELDGRATING_H
#define FULLFIELDGRATING_H

#include "grating.h"



class FullFieldGrating : public Grating
{
public:
    FullFieldGrating(Integrator *integrator, vec2 kd, double wd, double contrast);
    ~FullFieldGrating();

    // Stimulus interface
private:
    virtual double valueAtPoint(vec2 rVec, double t);
    virtual double fourierTransformAtFrequency(vec2 k, double w);
};

FullFieldGrating createFullFieldGratingStimulus(Integrator *integrator, const Config *cfg);

#endif // FULLFIELDGRATING_H
