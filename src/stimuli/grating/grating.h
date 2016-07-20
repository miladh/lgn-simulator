#ifndef GRATING_H
#define GRATING_H

#include "stimuli/stimulus.h"
#include <memory>

#include <cstddef>
#include <memory>
#include <type_traits>
#include <utility>

namespace lgnSimulator {
class Grating : public Stimulus
{
public:
    Grating(Integrator * const integrator,
            double spatialFreq, double orientation, double temporalFreq,
            double contrast, double maskSize = 0.0, double phase=0.0);
    ~Grating();

    // Stimulus interface
public:
    virtual void computeSpatiotemporal();
    virtual void computeFourierTransform();



    double spatialFreq() const;
    double orientation(bool inDegrees=false) const;
    double contrast() const;
    double maskSize() const;
    double temporalFreq() const;
    vec2 kVec() const;
    string mask() const;



protected:
    vec2 m_kVec = {0,0};
    double m_k = 0.0;
    double m_orientation = 0.0;
    double m_w = 0;
    double m_contrast = 0.0;
    double m_maskSize = 0.0;
    double m_phase = 0.0;
    string m_mask;

    virtual double valueAtPoint(vec2 r, double t) const  = 0;
    virtual complex<double> fourierTransformAtFrequency(vec2 k, double w) const= 0;

private:
    void setSpatialFreq(double spatialFreq);
    void setOrientation(double orientation);
};

}
std::unique_ptr<lgnSimulator::Grating> createGratingStimulus(
        lgnSimulator::Integrator* const integrator,
        const YAML::Node &cfg);

#endif // GRATING_H
