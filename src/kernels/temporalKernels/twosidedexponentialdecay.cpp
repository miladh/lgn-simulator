#include "twosidedexponentialdecay.h"

/*!
  \class lgnSimulator::TwoSidedExponentialDecay
  \inmodule lgnSimulator
  \ingroup lgnSimulator-temporalKernel
  \brief Temporal two sided exponential decay kernel.
 */


using namespace lgnSimulator;

TwoSidedExponentialDecay::TwoSidedExponentialDecay(double tau, double delay)
    : m_tau(tau)
    , m_delay(delay)
{

}


double lgnSimulator::TwoSidedExponentialDecay::temporal(double t) const
{
    return exp(-fabs(t - m_delay)/m_tau)/m_tau;
}

complex<double> lgnSimulator::TwoSidedExponentialDecay::fourierTransform(double w) const
{
    return 2. * exp(core::i*w*m_delay)/(1.0 + w*w * m_tau*m_tau);
}


TwoSidedExponentialDecay
createTemporalTwoSidedExponentialDecayKernel(const YAML::Node &cfg)
{

    double tau = cfg["tau"].as<double>();
    double delay = cfg["delay"].as<double>();

    return TwoSidedExponentialDecay(tau, delay);

}
