#include "spatiallyconstant.h"

using namespace lgnSimulator;

SpatiallyConstant::SpatiallyConstant(double constant, double spatialFreqResolution)
    : m_constant(constant)
    , m_peak(1./spatialFreqResolution/spatialFreqResolution)
{

}

double lgnSimulator::SpatiallyConstant::spatial(vec2 r) const
{
    (void)r;
    return m_constant;
}

complex<double> lgnSimulator::SpatiallyConstant::fourierTransform(vec2 k) const
{
    return 4.*core::pi*core::pi * m_constant * m_peak * Special::delta(k, vec2{0,0});
}


lgnSimulator::SpatiallyConstant createSpatiallyConstantKernel(const YAML::Node &cfg)
{
    double constant = cfg["constant"].as<double>();
    double ds = cfg["ds"].as<double>();
    double ns = cfg["ns"].as<double>();
    double peak = 2.*core::pi/pow(2,ns)/ds;

    return SpatiallyConstant(constant, peak);

}
