#include "dampedoscillator.h"

using namespace lgnSimulator;


DampedOscillator::DampedOscillator(double phaseDuration, double weight, double delay)
    : m_phaseDuration(phaseDuration)
    , m_weight(weight)
    , m_delay(delay)
{

}

DampedOscillator::~DampedOscillator()
{

}

double DampedOscillator::temporal(double t) const
{
    double tt = t - m_delay;
    if(tt <= m_phaseDuration){
        return sin(core::pi/m_phaseDuration * tt) * Special::heaviside(tt);
    }else if(tt <= 2*m_phaseDuration){
        return m_weight * sin(core::pi/m_phaseDuration * tt);
    }else{
        return 0.0;
    }
}

complex<double> DampedOscillator::fourierTransform(double w) const
{
    double factor = core::pi * m_phaseDuration/
            (core::pi*core::pi - m_phaseDuration * m_phaseDuration * w * w);
    complex<double> expTerm = exp(core::i * m_phaseDuration * w);
    complex<double> term1 = 1. + (1. - m_weight) * expTerm;
    complex<double> term2 = m_weight * exp(core::i * m_phaseDuration * 2.*w);

    return factor * exp(core::i * m_delay * w) * (term1 - term2);
}



DampedOscillator createTemporalDampedOscillatorKernel(const YAML::Node &cfg)
{
    double phaseDuration = cfg["phaseDuration"].as<double>();
    double weight = cfg["weight"].as<double>();
    double delay = cfg["delay"].as<double>();

    return DampedOscillator(phaseDuration, weight, delay);

}
