#include "decayingexponential.h"

DecayingExponential::DecayingExponential(double tau, double delay)
    : m_tau(tau)
    , m_delay(delay)
{

}

DecayingExponential::~DecayingExponential()
{

}

double DecayingExponential::real(double t)
{

    return 1./ m_tau * exp(-(t - m_delay)/m_tau) * Functions::heaviside(t - m_delay);
}

double DecayingExponential::complex(double w)
{
    return cos(w*m_delay)/ (1 + w*w * m_tau*m_tau);
}

