#include "decayingexponential.h"

using namespace lgnSimulator;


DecayingExponential::DecayingExponential(double tau, double delay)
    : m_tau(tau)
    , m_delay(delay)
{

}

DecayingExponential::~DecayingExponential()
{

}

double DecayingExponential::temporal(double t)
{

    return 1./ m_tau * exp(-(t - m_delay)/m_tau)
            * Functions::heaviside(t - m_delay);
}

complex<double> DecayingExponential::fourierTransform(double w)
{
//    return cos(w*m_delay)/ (1 + w*w * m_tau*m_tau);
        return exp(-complex<double>(0,1)*w*m_delay)
                /(complex<double>(1,0) + complex<double>(0,1)*w* m_tau);
}


DecayingExponential createTemporalDecayingExponentialKernel(const YAML::Node &cfg)
{

    double tau = cfg["tau"].as<double>();
    double delay = cfg["delay"].as<double>();

    return DecayingExponential(tau, delay);

}
