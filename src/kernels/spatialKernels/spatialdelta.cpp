#include "spatialdelta.h"

using namespace lgnSimulator;
SpatialDelta::SpatialDelta(double weight, double spatialResolution, vec2 shift)
    : m_weight(weight)
    , m_peak(1./spatialResolution)
    , m_shift(shift)
{
}


double lgnSimulator::SpatialDelta::spatial(vec2 r) const
{
    return m_weight * m_peak * m_peak
            * Special::delta(m_shift(0), r(0))
            * Special::delta(m_shift(1), r(1));
}

complex<double> lgnSimulator::SpatialDelta::fourierTransform(vec2 k) const
{
    return m_weight * exp(-core::i * dot(k, m_shift));
}



SpatialDelta createSpatialDeltaKernel(const YAML::Node &cfg)
{

    double weight = cfg["weight"].as<double>();
    double ds = cfg["ds"].as<double>();
    vec2 shift = {cfg["shift"][0].as<double>(),
                  cfg["shift"][1].as<double>()};


    return SpatialDelta(weight, ds, shift);
}
