#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <iostream>
#include <armadillo>
#include <math.h>
#include <fftw3.h>
#include <yaml-cpp/yaml.h>

#include "math/functions.h"
#include "math/ffthelper.h"

using namespace std;
using namespace arma;

namespace lgnSimulator {

class Integrator
{
public:
    Integrator(int nt, double dt, int ns, double ds);
    ~Integrator();

    cx_cube backwardFFT(cx_cube data);
    cx_mat backwardFFT(cx_mat data);

    cx_cube forwardFFT(cx_cube data);
    cx_mat forwardFFT(cx_mat data);

    vec timeVec() const;
    vec coordinateVec() const;
    vec temporalFreqVec() const;
    vec spatialFreqVec() const;

    int nPointsTemporal() const;
    int nPointsSpatial() const;
    double temporalFreqResolution() const;
    double spatialFreqResolution() const;
    double timeInterval() const;
    double lengthInterval() const;
    double temporalSamplingFreq() const;
    double spatialSamplingFreq() const;
    double dt() const;
    double ds() const;

private:
    int m_nPointsTemporal= 0;
    int m_nPointsSpatial = 0;

    double m_dt = 0;
    double m_ds = 0;

    double m_timeInterval = 0;
    double m_lengthInterval = 0;

    double m_temporalSamplingFreq = 0;
    double m_spatialSamplingFreq = 0;

    double m_dw = 0;
    double m_dk = 0;

    vec m_timeVec;
    vec m_coordinateVec;
    vec m_temporalFreqs;
    vec m_spatialFreqs;

};
}

lgnSimulator::Integrator createIntegrator(const YAML::Node *cfg);

#endif // INTEGRATOR_H
