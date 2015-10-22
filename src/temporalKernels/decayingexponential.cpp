#include "decayingexponential.h"

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

double DecayingExponential::fourierTransform(double w)
{
    return cos(w*m_delay)/ (1 + w*w * m_tau*m_tau);
}


DecayingExponential createDecayingExponentialTemporalKernel(const Config *cfg)
{
    const Setting & root = cfg->getRoot();
    double tau = root["temporalKernelSettings"]["tau"];
    double delay = root["temporalKernelSettings"]["delay"];

    return DecayingExponential(tau, delay);

}
