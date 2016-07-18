#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <iostream>
#include <armadillo>
#include <math.h>
#include <fftw3.h>
#include <yaml-cpp/yaml.h>

#include "helper/special.h"
#include "helper/ffthelper.h"
#include "helper/helperconstants.h"

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

    cube backwardFFT(cx_cube in);

    cx_cube forwardFFT(cube data);
//    cx_mat forwardFFT(mat data) const;

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
    const int m_nPointsTemporal= 0;
    const int m_nPointsSpatial = 0;

    const double m_temporalResolution = 0;
    const double m_spatialResolution = 0;

    const double m_timeInterval = 0;
    const double m_lengthInterval = 0;

    const double m_temporalSamplingFreq = 0;
    const double m_spatialSamplingFreq = 0;

    const double m_temporalFreqResolution = 0;
    const double m_spatialFreqResolution = 0;

    vec m_timeVec;
    vec m_spatialVec;
    vec m_temporalFreqs;
    vec m_spatialFreqs;

    cx_cube m_real;
    cx_cube m_complex;

    fftw_plan m_forwardPlan;
    fftw_plan m_backwardPlan;

};
}

lgnSimulator::Integrator createIntegrator(const YAML::Node &cfg);

#endif // INTEGRATOR_H
