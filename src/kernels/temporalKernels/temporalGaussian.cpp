#include "temporalGaussian.h"

using namespace lgnSimulator;

TemporalGaussian::TemporalGaussian(double A, double a, double delay)
    : m_A(A)
    , m_a(a)
    , m_delay(delay)
{

}

double TemporalGaussian::temporal(double t)
{
     return m_A / (m_a*m_a) / PI * exp(-(t-m_delay) / (m_a*m_a));
}

complex<double> TemporalGaussian::fourierTransform(double w)
{
    return m_A * exp(-w*w * m_a*m_a / 4.) * exp(complex<double>(0,1) * m_delay * w) ;
}


TemporalGaussian createTemporalGaussianKernel(const YAML::Node *cfg)
{
    double A = (*cfg)["A"].as<double>();
    double a = (*cfg)["a"].as<double>();
    double delay = (*cfg)["delay"].as<double>();

    return TemporalGaussian(A, a, delay);
}
