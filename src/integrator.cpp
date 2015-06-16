#include "integrator.h"

Integrator::Integrator(double upperLim, double lowerLim, double panels)
    : m_upperLim(upperLim)
    , m_lowerLim(lowerLim)
    , m_panels(panels)
    , m_stepLength((m_upperLim - m_lowerLim)/m_panels)

{

}

Integrator::~Integrator()
{

}
double Integrator::upperLim() const
{
    return m_upperLim;
}

double Integrator::lowerLim() const
{
    return m_lowerLim;
}

double Integrator::panels() const
{
    return m_panels;
}

double Integrator::stepLength() const
{
    return m_stepLength;
}




