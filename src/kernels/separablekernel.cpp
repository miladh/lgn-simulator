#include "separablekernel.h"

using namespace lgnSimulator;

SeparableKernel::SeparableKernel(SpatialKernel *spatialKernel,
                                 TemporalKernel *temporalKernel)
    : m_spatialKernel(spatialKernel)
    , m_temporalKernel(temporalKernel)

{

}

double SeparableKernel::spatiotemporal(vec2 r, double t)
{
    return m_spatialKernel->spatial(r) * m_temporalKernel->temporal(t);
}

complex<double> SeparableKernel::fourierTransform(vec2 k, double w)
{
    return m_spatialKernel->fourierTransform(k) * m_temporalKernel->fourierTransform(w);
}
