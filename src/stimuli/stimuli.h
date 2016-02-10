#ifndef STIMULI_H
#define STIMULI_H

#include <math.h>
#include <armadillo>
#include <yaml-cpp/yaml.h>

#include "integrator.h"

using namespace std;
using namespace arma;
namespace lgnSimulator {
class Stimulus
{
public:
    Stimulus(Integrator *integrator);
    ~Stimulus();

    virtual void computeSpatiotemporal() = 0;
    virtual void computeFourierTransform() = 0;

    cube spatioTemporal() const;
    cx_cube fourierTransform() const;

    void clearSpatioTemporal();
    void clearFourierTransform();

    string type() const;

protected:
    cube m_spatioTemporal;
    cx_cube m_fourierTransform;

    vec m_spatialVec;
    vec m_spatialFreqs;

    vec m_timeVec;
    vec m_temporalFreqs;

    Integrator *m_integrator;

    string m_type;


};
}
#endif // STIMULI_H
