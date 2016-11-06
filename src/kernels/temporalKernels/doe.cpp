#include "doe.h"

/*!
  \class lgnSimulator::DOE
  \inmodule lgnSimulator
  \ingroup lgnSimulator-temporalKernel
  \brief Temporal difference of exponentials kernel.
 */

using namespace lgnSimulator;

DOE::DOE(double cenLatency, double surLatency, double delay)
    : m_cenLatency(cenLatency)
    , m_surLatency(surLatency)
    , m_delay(delay)
{

}


double DOE::temporal(double t) const
{
    if(Special::heaviside(t - m_delay) == 0) {
        return 0;
    }else{
        double center = exp(-(t - m_delay)/m_cenLatency)
                /m_cenLatency/m_cenLatency * (t - m_delay);
        double surround = exp(-(t - m_delay)/m_surLatency)
                /m_surLatency/m_surLatency * (t - m_delay);
        return center - surround;
    }
}

complex<double> DOE::fourierTransform(double w) const
{
    complex<double> expFactor = exp(core::i*w*m_delay);
    complex<double> center = expFactor
                    /(1. - m_cenLatency*w*core::i)
                    /(1. - m_cenLatency*w*core::i);
    complex<double> surround = expFactor
                    /(1. - m_surLatency*w*core::i)
                    /(1. - m_surLatency*w*core::i);

    return center - surround;

}

DOE createTemporalDOEKernel(const YAML::Node &cfg)
{

    double cenLatency = cfg["cenLatency"].as<double>();
    double surLatency = cfg["surLatency"].as<double>();
    double delay = cfg["delay"].as<double>();

    return DOE(cenLatency, surLatency, delay);

}
