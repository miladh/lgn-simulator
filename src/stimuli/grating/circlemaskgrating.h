#ifndef CIRCLEMASKGRATING_H
#define CIRCLEMASKGRATING_H

#include "grating.h"


class CircleMaskGrating : public Grating
{
public:
    CircleMaskGrating(Integrator *integrator,
                 vec2 kd, double wd, double contrast, double maskSize);
    ~CircleMaskGrating();

private:
    virtual double valueAtPoint(vec2 rVec, double t);
    virtual double fourierTransformAtFrequency(vec2 kVec, double w);

};

CircleMaskGrating createCircleMaskGratingStimulus(Integrator *integrator, const Config* cfg);

#endif // CIRCLEMASKGRATING_H
