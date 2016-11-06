#include "biphasic.h"

/*!
  \class lgnSimulator::Biphasic
  \inmodule lgnSimulator
  \ingroup lgnSimulator-temporalKernel
  \brief Temporal biphasic kernel.
 */

using namespace lgnSimulator;


Biphasic::Biphasic(double phaseDuration, double dampingFactor, double delay)
    : m_phaseDuration(phaseDuration)
    , m_dampingFactor(dampingFactor)
    , m_delay(delay)
{

}

Biphasic::~Biphasic()
{

}

double Biphasic::temporal(double t) const
{
    double tt = t - m_delay;
    if(Special::heaviside(tt) == 0){
        return 0;
    }

    if(tt <= m_phaseDuration){
        return sin(core::pi/m_phaseDuration * tt);
    }else if(tt <= 2*m_phaseDuration){
        return m_dampingFactor * sin(core::pi/m_phaseDuration * tt);
    }else{
        return 0.0;
    }
}

complex<double> Biphasic::fourierTransform(double w) const
{
    double factor = core::pi * m_phaseDuration/
            (core::pi*core::pi - m_phaseDuration * m_phaseDuration * w * w);
    complex<double> expTerm = exp(core::i * m_phaseDuration * w);
    complex<double> term1 = 1. + (1. - m_dampingFactor) * expTerm;
    complex<double> term2 = m_dampingFactor * exp(core::i * m_phaseDuration * 2.*w);

    return factor * exp(core::i * m_delay * w) * (term1 - term2);
}



Biphasic createTemporalBiphasicKernel(const YAML::Node &cfg)
{
    double phaseDuration = cfg["phaseDuration"].as<double>();
    double weight = cfg["weight"].as<double>();
    double delay = cfg["delay"].as<double>();

    return Biphasic(phaseDuration, weight, delay);

}
