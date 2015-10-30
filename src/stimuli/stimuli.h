#ifndef STIMULI_H
#define STIMULI_H

#include <math.h>
#include <armadillo>
#include <libconfig.h++>

#include "integrator.h"

using namespace std;
using namespace arma;
using namespace libconfig;

class Stimulus
{
public:
    Stimulus(Integrator *integrator);
    ~Stimulus();

    virtual void computeSpatiotemporal() = 0;
    virtual void computeFourierTransform() = 0;

    cube spatioTemporal() const;
    cx_cube fourierTransform() const;

protected:
    cube m_spatioTemporal;
    cx_cube m_fourierTransform;

    vec m_coordinateVec;
    vec m_spatialFreqs;

    vec m_timeVec;
    vec m_temporalFreqs;

    Integrator *m_integrator;

    void computeSpatiotemporalAnalytic();
    void computeFourierTransformAnalytic();

private:
    virtual double valueAtPoint(vec2 rVec, double t) = 0;
    virtual double fourierTransformAtFrequency(vec2 k, double w) = 0;

};

#endif // STIMULI_H
