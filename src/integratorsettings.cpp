#include "integratorsettings.h"

IntegratorSettings::IntegratorSettings(int nt, int ns, double maxTime)
    : m_nPointsTemporal(pow(2,nt))
    , m_nPointsSpatial(pow(2,ns))
    , m_maxTime(maxTime)
{

}

IntegratorSettings::~IntegratorSettings()
{

}
int IntegratorSettings::nPointsTemporal() const
{
    return m_nPointsTemporal;
}

int IntegratorSettings::nPointsSpatial() const
{
    return m_nPointsSpatial;
}

double IntegratorSettings::maxTime() const
{
    return m_maxTime;
}





