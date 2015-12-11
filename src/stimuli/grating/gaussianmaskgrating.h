#ifndef GAUSSIANMASKGRATING_H
#define GAUSSIANMASKGRATING_H

#include "grating.h"

class GaussianMaskGrating : public Grating
{
public:
    GaussianMaskGrating(Integrator *integrator,
                        vec2 kd, double wd, double contrast, double maskSize);
    ~GaussianMaskGrating();

private:
    double valueAtPoint(vec2 rVec, double t);
    double fourierTransformAtFrequency(vec2 kVec, double w);

};

GaussianMaskGrating createGaussianMaskGratingStimulus(Integrator *integrator, const Config* cfg);

#endif // GAUSSIANMASKGRATING_H
