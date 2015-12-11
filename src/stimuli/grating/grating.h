#ifndef GRATING_H
#define GRATING_H

#include "../stimuli.h"


class Grating : public Stimulus
{
public:
    Grating(Integrator *integrator,
            vec2 kd, double wd, double contrast,
            string mask = "None", double maskSize = 0.0);
    ~Grating();

    // Stimulus interface
public:
    virtual void computeSpatiotemporal();
    virtual void computeFourierTransform();

private:
    virtual double valueAtPoint(vec2 rVec, double t) = 0;
    virtual double fourierTransformAtFrequency(vec2 kVec, double w) = 0;

protected:
    vec2 m_k = {0,0};
    double m_w = 0;
    double m_contrast = 0.0;
    double m_maskSize = 0.0;

};

Grating* setGratingType(Integrator *integrator,
                        vec2 kd, double wd, double contrast,
                        string mask, double maskSize);

Grating *createGratingStimulus(Integrator *integrator, const Config* cfg);

#endif // GRATING_H
