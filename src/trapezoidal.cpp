#include "trapezoidal.h"

Trapezoidal::Trapezoidal(double upperLim, double lowerLim, double steps)
    : Integrator(upperLim, lowerLim, steps)
{

}

Trapezoidal::~Trapezoidal()
{

}


void Trapezoidal::integrate()
{

    m_result = 0.;
    for (int i = 0; i < panels(); i++){
        m_result += advance(lowerLim() + stepLength()*i);
    }

    m_result *=stepLength();
}


double Trapezoidal::advance(double x) const
{
    return (x*x + pow(x+stepLength(),2)) * 0.5;
}


double Trapezoidal::result() const
{
    return m_result;
}




