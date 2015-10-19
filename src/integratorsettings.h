#ifndef INTEGRATORSETTINGS_H
#define INTEGRATORSETTINGS_H

#include <armadillo>
#include <iostream>

using namespace std;
using namespace arma;

class IntegratorSettings
{
public:
    IntegratorSettings(int nt, double dt, int ns, double ds);
    ~IntegratorSettings();

    int nPointsTemporal() const;
    int nPointsSpatial() const;
    double temporalResolution() const;
    double spatialResolution() const;

private:
    int m_nPointsTemporal = 0;
    int m_nPointsSpatial = 0;
    double m_dt = 0;
    double m_ds = 0;
};

#endif // INTEGRATORSETTINGS_H
