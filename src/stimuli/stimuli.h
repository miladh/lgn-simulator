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

    void computeSpatiotemporal();
    void computeFourierTransform();
    void setSpatial(cube valueAtPoint);

    cube spatioTemporal() const;
    cx_cube fourierTransform() const;

protected:
    cube m_spatioTemporal;
    cx_cube m_fourierTransform;

    vec m_coordinateVec;
    vec m_spatialFreqs;

    vec timeVec;
    vec m_temporalFreqs;

    Integrator *m_integrator;

private:
    virtual double valueAtPoint(vec2 rVec, double t) = 0;
    virtual double fourierTransformAtFrequency(vec2 k, double w) = 0;

};

#endif // STIMULI_H
