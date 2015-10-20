#include "integratorsettings.h"

IntegratorSettings::IntegratorSettings(int nt, double dt, int ns, double ds)
    : m_nPointsTemporal(pow(2,nt))
    , m_dt(dt)
    , m_nPointsSpatial(pow(2,ns))
    , m_ds(ds)
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

double IntegratorSettings::temporalResolution() const
{
    return m_dt;
}

double IntegratorSettings::spatialResolution() const
{
    return m_ds;
}








