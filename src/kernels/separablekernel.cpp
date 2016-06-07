#include "separablekernel.h"

using namespace lgnSimulator;

SeparableKernel::SeparableKernel(double weight,
                                 SpatialKernel *spatialKernel,
                                 TemporalKernel *temporalKernel)
    : Kernel(weight)
    , m_spatialKernel(spatialKernel)
    , m_temporalKernel(temporalKernel)

{

}

double SeparableKernel::spatiotemporal(vec2 r, double t) const
{
    return m_weight * m_spatialKernel->spatial(r) * m_temporalKernel->temporal(t);
}

complex<double> SeparableKernel::fourierTransform(vec2 k, double w) const
{
    return m_weight * m_spatialKernel->fourierTransform(k) * m_temporalKernel->fourierTransform(w);
}
