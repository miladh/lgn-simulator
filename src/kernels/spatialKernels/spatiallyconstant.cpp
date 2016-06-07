#include "spatiallyconstant.h"

using namespace lgnSimulator;


SpatiallyConstant::SpatiallyConstant(double lengthInterval, double spatialFreqResolution)
    : m_lengthInterval(lengthInterval)
    , m_peak(1./spatialFreqResolution/spatialFreqResolution)
{

}

double lgnSimulator::SpatiallyConstant::spatial(vec2 r) const
{
    (void)r;
    return 1./m_lengthInterval/m_lengthInterval;
}

complex<double> lgnSimulator::SpatiallyConstant::fourierTransform(vec2 k) const
{
    return 4.*core::pi*core::pi * m_peak / m_lengthInterval / m_lengthInterval
            * Special::delta(k, vec2{0,0});
}


lgnSimulator::SpatiallyConstant createSpatiallyConstantKernel(const YAML::Node &cfg)
{
    double ds = cfg["ds"].as<double>();
    double ns = cfg["ns"].as<double>();

    double lengthInterval =pow(2,ns)*ds;
    double peak = 2.*core::pi/lengthInterval;

    return SpatiallyConstant(lengthInterval, peak);

}
