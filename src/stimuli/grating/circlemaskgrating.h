#ifndef CIRCLEMASKGRATING_H
#define CIRCLEMASKGRATING_H

#include "grating.h"


class CircleMaskGrating : public Grating
{
public:
    CircleMaskGrating(Integrator *integrator,
                 vec2 kd, double wd, double contrast, double maskSize);
    ~CircleMaskGrating();

        // Stimulus interface
private:
    virtual double valueAtPoint(vec2 rVec, double t);
    virtual double fourierTransformAtFrequency(vec2 kVec, double w);

};


#endif // CIRCLEMASKGRATING_H
