#include "spatialgaussian.h"


using namespace lgnSimulator;

SpatialGaussian::SpatialGaussian(double a)
        : m_a(a)
{

}


double SpatialGaussian::spatial(vec2 r) const
{
    return 1. / (m_a*m_a) / core::pi * exp(-dot(r,r) / (m_a*m_a));
}

complex<double> SpatialGaussian::fourierTransform(vec2 k) const
{

    return exp(-dot(k,k) * m_a*m_a / 4.);
}

SpatialGaussian createSpatialGaussianKernel(const YAML::Node &cfg)
{
    double a = cfg["a"].as<double>();

    return SpatialGaussian(a);
}
