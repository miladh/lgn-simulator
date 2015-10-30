#ifndef OSCILLATINGGAUSSIAN_H
#define OSCILLATINGGAUSSIAN_H

#include "stimuli.h"



class OscillatingGaussian : public Stimulus
{
public:
    OscillatingGaussian(Integrator *integrator,
                        double exponent, double wd);
    ~OscillatingGaussian();

    // Stimulus interface
public:
    virtual void computeSpatiotemporal();
    virtual void computeFourierTransform();
private:
    double valueAtPoint(vec2 rVec, double t);
    double fourierTransformAtFrequency(vec2 kVec, double w);


private:
    double m_exponent;
    double m_wd;

};


OscillatingGaussian createOscillatingGaussianStimulus(Integrator *integrator,
                                                      const Config *cfg);

#endif // OSCILLATINGGAUSSIAN_H
