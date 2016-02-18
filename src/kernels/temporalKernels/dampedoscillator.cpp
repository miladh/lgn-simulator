#include "dampedoscillator.h"

using namespace lgnSimulator;


DampedOscillator::DampedOscillator(double phaseDuration, double weight)
    : m_phaseDuration(phaseDuration)
    , m_weight(weight)
{

}

DampedOscillator::~DampedOscillator()
{

}

double DampedOscillator::temporal(double t) const
{
    if(t>= 0 && t<=m_phaseDuration){
        return sin(core::pi/m_phaseDuration * t);
    }else if(t > m_phaseDuration && t <= 2*m_phaseDuration){
        return m_weight * sin(core::pi/m_phaseDuration * t);
    }else{
        return 0.0;
    }
}

complex<double> DampedOscillator::fourierTransform(double w) const
{
    double factor = core::pi*m_phaseDuration/
            (core::pi*core::pi - m_phaseDuration * m_phaseDuration * w * w);


    return -factor*((complex<double>(1,0)+ exp(core::i*m_phaseDuration*w))
            * (complex<double>(-1,0)
               + m_weight *exp(core::i*m_phaseDuration*w)));
}



DampedOscillator createTemporalDampedOscillatorKernel(const YAML::Node &cfg)
{
    double phaseDuration = cfg["phaseDuration"].as<double>();
    double weight = cfg["weight"].as<double>();

    return DampedOscillator(phaseDuration, weight);

}
