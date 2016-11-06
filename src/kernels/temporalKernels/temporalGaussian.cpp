#include "temporalGaussian.h"

/*!
  \class lgnSimulator::TemporalGaussian
  \inmodule lgnSimulator
  \ingroup lgnSimulator-temporalKernel
  \brief Temporal Gaussian kernel.
 */

using namespace lgnSimulator;

TemporalGaussian::TemporalGaussian(double A, double a, double delay)
    : m_A(A)
    , m_a(a)
    , m_delay(delay)
{

}

double TemporalGaussian::temporal(double t) const
{
     return m_A / (m_a*sqrt(core::pi)) * exp(-(t-m_delay)*(t-m_delay) / (m_a*m_a));
}

complex<double> TemporalGaussian::fourierTransform(double w) const
{
    return m_A * exp(-w*w * m_a*m_a / 4.) * exp(core::i * m_delay * w) ;
}


TemporalGaussian createTemporalGaussianKernel(const YAML::Node &cfg)
{
    double A = cfg["A"].as<double>();
    double a = cfg["a"].as<double>();
    double delay = cfg["delay"].as<double>();

    return TemporalGaussian(A, a, delay);
}
