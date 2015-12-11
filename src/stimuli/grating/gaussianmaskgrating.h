#ifndef GAUSSIANMASKGRATING_H
#define GAUSSIANMASKGRATING_H

#include "grating.h"
#include "../../spatialKernels/dog.h"


class GaussianMaskGrating : public Grating
{
public:
    GaussianMaskGrating(Integrator *integrator,
                        vec2 kd, double wd, double contrast, double maskSize);
    ~GaussianMaskGrating();

    // Stimulus interface
private:
    virtual double valueAtPoint(vec2 rVec, double t);
    virtual double fourierTransformAtFrequency(vec2 kVec, double w);

    DOG* m_dog;

};

GaussianMaskGrating createGaussianMaskGratingStimulus(Integrator *integrator, const Config* cfg);

#endif // GAUSSIANMASKGRATING_H
