#include "gaussian.h"


using namespace lgnSimulator;

Gaussian::Gaussian(double A, double a)
        : m_A(A)
        , m_a(a)
{

}


double Gaussian::spatial(vec2 r)
{
    return m_A / (m_a*m_a) / PI * exp(-dot(r,r) / (m_a*m_a));
}

complex<double> Gaussian::fourierTransform(vec2 k)
{

    return m_A * exp(-dot(k,k) * m_a*m_a / 4.);
}



Gaussian createSpatialGaussianKernel(const YAML::Node *cfg)
{
    double A = (*cfg)["A"].as<double>();
    double a = (*cfg)["a"].as<double>();


    return Gaussian(A, a);
}
