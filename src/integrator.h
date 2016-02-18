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
    Integrator(const int nt,
               const double temporalResolution,
               const int ns,
               const double spatialResolution);
    ~Integrator();

    cx_cube backwardFFT(cx_cube data) const;
    cx_mat backwardFFT(cx_mat data) const;

    cx_cube forwardFFT(cx_cube data) const;
    cx_mat forwardFFT(cx_mat data) const;

    vec timeVec() const;
    vec spatialVec() const;
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
    double temporalResolution() const;
    double spatialResolution() const;

private:
    int m_nPointsTemporal= 0;
    int m_nPointsSpatial = 0;

    double m_temporalResolution = 0;
    double m_spatialResolution = 0;

    double m_timeInterval = 0;
    double m_lengthInterval = 0;

    double m_temporalSamplingFreq = 0;
    double m_spatialSamplingFreq = 0;

    double m_temporalFreqResolution = 0;
    double m_spatialFreqResolution = 0;

    vec m_timeVec;
    vec m_spatialVec;
    vec m_temporalFreqs;
    vec m_spatialFreqs;

};
}

lgnSimulator::Integrator createIntegrator(const YAML::Node &cfg);

#endif // INTEGRATOR_H
