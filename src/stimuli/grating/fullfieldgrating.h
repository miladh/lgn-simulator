#ifndef FULLFIELDGRATING_H
#define FULLFIELDGRATING_H

#include "grating.h"


namespace lgnSimulator {
class FullFieldGrating : public Grating
{
public:
    FullFieldGrating(const Integrator &integrator, vec2 kd, double wd, double contrast);
    ~FullFieldGrating();

    // Grating interface
private:
    virtual double valueAtPoint(vec2 rVec, double t);
    virtual complex<double> fourierTransformAtFrequency(vec2 k, double w);
};

}
#endif // FULLFIELDGRATING_H
