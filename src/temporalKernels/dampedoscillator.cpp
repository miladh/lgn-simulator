#include "dampedoscillator.h"


DampedOscillator::DampedOscillator(double phaseDuration, double weight)
    : m_phaseDuration(phaseDuration)
    , m_weight(weight)
{

}

DampedOscillator::~DampedOscillator()
{

}

double DampedOscillator::real(double t)
{
    if(t>= 0 || t<=m_phaseDuration){
        return sin(PI/m_phaseDuration * t);
    }else if(t > m_phaseDuration || t <= 2*m_phaseDuration){
        return m_weight * sin(PI/m_phaseDuration * t);
    }else{
        return 0;
    }
}

double DampedOscillator::complex(double w)
{
    double factor = PI*m_phaseDuration/
            (PI*PI - m_phaseDuration * m_phaseDuration * w * w);

    return factor * (1 + cos(m_phaseDuration * w)*(1 - m_weight)
                     - m_weight * cos(2* m_phaseDuration * w));
}
