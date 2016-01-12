#ifndef GRATING_H
#define GRATING_H

#include "../stimuli.h"
#include <memory>

#include <cstddef>
#include <memory>
#include <type_traits>
#include <utility>

class Grating : public Stimulus
{
public:
    Grating(Integrator *integrator,
            vec2 kd, double wd, double contrast, double maskSize = 0.0);
    ~Grating();

    // Stimulus interface
public:
    virtual void computeSpatiotemporal();
    virtual void computeFourierTransform();


protected:
    vec2 m_k = {0,0};
    double m_w = 0;
    double m_contrast = 0.0;
    double m_maskSize = 0.0;

    virtual double valueAtPoint(vec2 rVec, double t) = 0;
    virtual double fourierTransformAtFrequency(vec2 k, double w) = 0;

};

std::unique_ptr<Grating> createGratingStimulus(Integrator *integrator, const Config* cfg);

#endif // GRATING_H
