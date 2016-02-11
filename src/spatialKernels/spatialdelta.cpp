#include "spatialdelta.h"

using namespace lgnSimulator;
SpatialDelta::SpatialDelta(double weight, vec2 r0)
    : m_weight(weight)
    , m_r0(r0)
{

}


double lgnSimulator::SpatialDelta::spatial(vec2 r)
{
    return m_weight * Functions::delta(m_r0(0), r(0)) * Functions::delta(m_r0(0), r(1));
}

complex<double> lgnSimulator::SpatialDelta::fourierTransform(vec2 k)
{
    complex<double> i = complex<double>(0,1);
    return m_weight * exp(-i * dot(k, m_r0));
}



SpatialDelta createSpatialDeltaKernel(const YAML::Node *cfg)
{

    double weight = (*cfg)["weight"].as<double>();
    vec2 r0 = {0,0};

    return SpatialDelta(weight, r0);
}
