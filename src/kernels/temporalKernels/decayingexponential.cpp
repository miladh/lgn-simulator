#include "decayingexponential.h"

/*!
  \class lgnSimulator::DecayingExponential
  \inmodule lgnSimulator
  \ingroup lgnSimulator-temporalKernel
  \brief Temporal decaying exponential kernel.
 */

using namespace lgnSimulator;


DecayingExponential::DecayingExponential(double tau, double delay)
    : m_tau(tau)
    , m_delay(delay)
{

}

DecayingExponential::~DecayingExponential()
{

}

double DecayingExponential::temporal(double t) const
{

    return exp(-(t - m_delay)/m_tau)
            * Special::heaviside(t - m_delay)/m_tau;
}

complex<double> DecayingExponential::fourierTransform(double w) const
{
    return exp(core::i*w*m_delay)/(complex<double>(1,0) - core::i*w* m_tau);
}


DecayingExponential createTemporalDecayingExponentialKernel(const YAML::Node &cfg)
{

    double tau = cfg["tau"].as<double>();
    double delay = cfg["delay"].as<double>();

    return DecayingExponential(tau, delay);

}
