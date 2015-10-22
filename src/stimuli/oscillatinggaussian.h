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
private:
    virtual double valueAtPoint(vec2 rVec, double t);
    virtual double fourierTransformAtFrequency(vec2 k, double w);

    double m_exponent;
    double m_wd;
};


OscillatingGaussian createOscillatingGaussianStimulus(Integrator *integrator,
                                                      const Config *cfg);

#endif // OSCILLATINGGAUSSIAN_H
