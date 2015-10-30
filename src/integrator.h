#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <iostream>
#include <armadillo>
#include <math.h>
#include <libconfig.h++>
#include <fftw3.h>

#include "math/functions.h"
#include "math/ffthelper.h"

using namespace std;
using namespace arma;
using namespace libconfig;

class Integrator
{
public:
    Integrator(int nt, double dt, int ns, double ds);
    ~Integrator();

    cx_cube backwardFFT(cx_cube data);
    cx_mat backwardFFT(cx_mat data);

    vec timeVec() const;
    vec coordinateVec() const;
    vec temporalFreqVec() const;
    vec spatialFreqVec() const;

    int nPointsTemporal() const;
    int nPointsSpatial() const;

    double dt() const;
    double ds() const;

    double temporalFreqResolution() const;
    double spatialFreqResolution() const;

private:
    int m_nPointsTemporal= 0;
    int m_nPointsSpatial = 0;
    double m_maxT = 0;

    double m_dt = 0;
    double m_dw = 0;

    double m_ds = 0;
    double m_dk = 0;

    vec m_timeVec;
    vec m_coordinateVec;
    vec m_temporalFreqs;
    vec m_spatialFreqs;

};

Integrator createIntegrator(const Config* cfg);

#endif // INTEGRATOR_H
