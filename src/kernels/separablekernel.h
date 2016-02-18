#ifndef SEPARABLEKERNEL_H
#define SEPARABLEKERNEL_H

#include "kernel.h"
#include "spatialKernels/spatialkernel.h"
#include "temporalKernels/temporalkernel.h"

namespace lgnSimulator {
class SeparableKernel : public Kernel
{
public:
    SeparableKernel(SpatialKernel *spatialKernel,
                    TemporalKernel *temporalKernel);

    // Kernel interface
    virtual double spatiotemporal(vec2 r, double t) const;
    virtual complex<double> fourierTransform(vec2 k, double w) const;

private:
    SpatialKernel *m_spatialKernel;
    TemporalKernel *m_temporalKernel;
};

}
#endif // SEPARABLEKERNEL_H
