#include "spatialgaussian.h"


using namespace lgnSimulator;

SpatialGaussian::SpatialGaussian(double A, double a)
        : m_A(A)
        , m_a(a)
{

}


double SpatialGaussian::spatial(vec2 r) const
{
    return m_A / (m_a*m_a) / core::pi * exp(-dot(r,r) / (m_a*m_a));
}

complex<double> SpatialGaussian::fourierTransform(vec2 k) const
{

    return m_A * exp(-dot(k,k) * m_a*m_a / 4.);
}

SpatialGaussian createSpatialGaussianKernel(const YAML::Node &cfg)
{
    double A = cfg["A"].as<double>();
    double a = cfg["a"].as<double>();

    return SpatialGaussian(A, a);
}
