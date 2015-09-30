#ifndef INTEGRATORSETTINGS_H
#define INTEGRATORSETTINGS_H

#include <armadillo>
#include <iostream>

using namespace std;
using namespace arma;

class IntegratorSettings
{
public:
    IntegratorSettings(int nt, int ns, double maxTime);
    ~IntegratorSettings();

    int nPointsTemporal() const;
    int nPointsSpatial() const;
    double maxTime() const;

private:
    int m_nPointsTemporal = 0;
    int m_nPointsSpatial = 0;
    double m_maxTime = 0;
};

#endif // INTEGRATORSETTINGS_H
