#ifndef GAUSSIANMASKGRATING_H
#define GAUSSIANMASKGRATING_H

#include "grating.h"
#include "kernels/spatialKernels/spatialgaussian.h"

namespace lgnSimulator {

class GaussianMaskGrating : public Grating
{
public:
    GaussianMaskGrating(Integrator* const integrator,
                        double spatialFreq, double temporalFreq,
                        double contrast, double phase,
                        double orientation, double maskSize);
    ~GaussianMaskGrating();

    // Grating interface
protected:
    virtual double valueAtPoint(vec2 r, double t) const override;
    virtual complex<double> fourierTransformAtFrequency(vec2 k, double w) const override;

private:
    SpatialGaussian *m_gaussianMask = nullptr;
};
}

#endif // GAUSSIANMASKGRATING_H
