#ifndef GAUSSIANMASKGRATING_H
#define GAUSSIANMASKGRATING_H

#include "grating.h"
#include "spatialKernels/dog.h"

namespace edog {
class GaussianMaskGrating : public Grating
{
public:
    GaussianMaskGrating(Integrator *integrator,
                        vec2 kd, double wd, double contrast, double maskSize);
    ~GaussianMaskGrating();

    // Grating interface
private:
    virtual double valueAtPoint(vec2 rVec, double t);
    virtual complex<double> fourierTransformAtFrequency(vec2 kVec, double w);

    DOG* m_dog;

};
}
#endif // GAUSSIANMASKGRATING_H
