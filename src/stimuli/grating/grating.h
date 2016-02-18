#ifndef GRATING_H
#define GRATING_H

#include "stimuli/stimuli.h"
#include <memory>

#include <cstddef>
#include <memory>
#include <type_traits>
#include <utility>
namespace lgnSimulator {
class Grating : public Stimulus
{
public:
    Grating(const Integrator &integrator,
            vec2 kd, double wd, double contrast, double maskSize = 0.0);
    ~Grating();

    // Stimulus interface
public:
    virtual void computeSpatiotemporal();
    virtual void computeFourierTransform();


    double contrast() const;
    double maskSize() const;
    double w() const;
    vec2 k() const;
    string mask() const;



protected:
    vec2 m_k = {0,0};
    double m_w = 0;
    double m_contrast = 0.0;
    double m_maskSize = 0.0;
    string m_mask;

    virtual double valueAtPoint(vec2 rVec, double t) const = 0;
    virtual complex<double> fourierTransformAtFrequency(vec2 k, double w) const= 0;

};

}
std::unique_ptr<lgnSimulator::Grating> createGratingStimulus(
        const lgnSimulator::Integrator &integrator,
        const YAML::Node &cfg);

#endif // GRATING_H
