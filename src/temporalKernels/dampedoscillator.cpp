#include "dampedoscillator.h"


DampedOscillator::DampedOscillator(double phaseDuration, double weight)
    : m_phaseDuration(phaseDuration)
    , m_weight(weight)
{

}

DampedOscillator::~DampedOscillator()
{

}

double DampedOscillator::temporal(double t)
{
    if(t>= 0 && t<=m_phaseDuration){
        return sin(PI/m_phaseDuration * t);
    }else if(t > m_phaseDuration && t <= 2*m_phaseDuration){
        return m_weight * sin(PI/m_phaseDuration * t);
    }else{
        return 0.0;
    }
}

complex<double> DampedOscillator::fourierTransform(double w)
{
    double factor = PI*m_phaseDuration/
            (PI*PI - m_phaseDuration * m_phaseDuration * w * w);


    return -factor*((complex<double>(1,0)+ exp(complex<double>(0,1)*m_phaseDuration*w))
            * (complex<double>(-1,0)
               + m_weight *exp(complex<double>(0,1)*m_phaseDuration*w)));
}



DampedOscillator createDampedOscillatorTemporalKernel(const Config *cfg)
{
    const Setting & root = cfg->getRoot();
    double phaseDuration = root["temporalKernelSettings"]["phaseDuration"];
    double weight = root["temporalKernelSettings"]["weight"];

    return DampedOscillator(phaseDuration, weight);

}
