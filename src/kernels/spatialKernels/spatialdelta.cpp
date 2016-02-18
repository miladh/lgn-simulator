#include "spatialdelta.h"

using namespace lgnSimulator;
SpatialDelta::SpatialDelta(double weight, vec2 r0)
    : m_weight(weight)
    , m_r0(r0)
{

}


double lgnSimulator::SpatialDelta::spatial(vec2 r) const
{
    return m_weight * SpecialFunctions::delta(m_r0(0), r(0)) * SpecialFunctions::delta(m_r0(0), r(1));
}

complex<double> lgnSimulator::SpatialDelta::fourierTransform(vec2 k) const
{
    return m_weight * exp(-core::i * dot(k, m_r0));
}



SpatialDelta createSpatialDeltaKernel(const YAML::Node &cfg)
{

    double weight = cfg["weight"].as<double>();
    vec2 r0 = {0,0};

    return SpatialDelta(weight, r0);
}
