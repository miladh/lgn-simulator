#ifndef PATCHGRATING_H
#define PATCHGRATING_H

#include "stimuli.h"


class PatchGrating : public Stimulus
{
public:
    PatchGrating(Integrator *integrator, vec2 kd,
                 double wd, double contrast, double spotDiameterRatio);
    ~PatchGrating();


    // Stimulus interface
public:
    virtual void computeSpatiotemporal();
    virtual void computeFourierTransform();
private:
    double valueAtPoint(vec2 rVec, double t);
    double fourierTransformAtFrequency(vec2 kVec, double w);

private:
    vec2 m_k = {0,0};
    double m_w = 0;
    double m_contrast = 0.0;
    double m_spotDiameter = 0.0;
};

PatchGrating createPatchGratingStimulus(Integrator *integrator, const Config* cfg);

#endif // PATCHGRATING_H
