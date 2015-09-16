#ifndef STIMULI_H
#define STIMULI_H

#include <math.h>
#include <libconfig.h++>
#include <armadillo>
#include <fftw3.h>
#include <boost/math/special_functions/bessel.hpp>

#include "../math/functions.h"

using namespace std;
using namespace libconfig;
using namespace arma;

class Stimuli
{
public:
    Stimuli(const Config *cfg);
    ~Stimuli();

    virtual double spatial(vec2 rVec, double t) = 0;
    virtual double frequency(vec2 k, double w) = 0;

    void computeSpatial();
    void computeFrequency();

    cube spatial() const;
    cx_cube frequency() const;

    double w() const;

    void setSpatial(cube spatial);

protected:
    int m_nPoints = 0;
    int m_nSteps = 0;
    double m_w = 0;
    cube m_spatial;
    cx_cube m_frequency;

    vec2 m_k = {0,0};

    vec m_spatialMesh;
    vec m_spatialFreqs;

    vec m_temporalMesh;
    vec m_temporalFreqs;


};

#endif // STIMULI_H
