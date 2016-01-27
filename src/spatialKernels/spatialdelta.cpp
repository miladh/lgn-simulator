#include "spatialdelta.h"

using namespace edog;
SpatialDelta::SpatialDelta(double weight, vec2 r0)
    : m_weight(weight)
    , m_r0(r0)
{

}


double edog::SpatialDelta::spatial(vec2 rVec)
{
    return m_weight * Functions::delta(m_r0(0), rVec(0)) * Functions::delta(m_r0(0), rVec(1));
}

complex<double> edog::SpatialDelta::fourierTransform(vec2 kVec)
{
    complex<double> i = complex<double>(0,1);
    return exp(-i * dot(kVec, m_r0));
}



SpatialDelta createSpatialDeltaSpatialKernel(const YAML::Node *cfg)
{

    double weight = (*cfg)["weight_delta"].as<double>();
    vec2 r0 = {0,0};

    return SpatialDelta(weight, r0);
}
