#ifndef CIRCLEMASKGRATING_H
#define CIRCLEMASKGRATING_H

#include "grating.h"

namespace edog {
class CircleMaskGrating : public Grating
{
public:
    CircleMaskGrating(Integrator *integrator,
                 vec2 kd, double wd, double contrast, double maskSize);
    ~CircleMaskGrating();

        // Grating interface
private:
    virtual double valueAtPoint(vec2 rVec, double t);
    virtual complex<double> fourierTransformAtFrequency(vec2 kVec, double w);

};

}
#endif // CIRCLEMASKGRATING_H
