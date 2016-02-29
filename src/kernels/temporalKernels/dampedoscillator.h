#ifndef DAMPEDOSCILLATOR_H
#define DAMPEDOSCILLATOR_H


#include "kernels/temporalKernels/temporalkernel.h"

namespace lgnSimulator {

class DampedOscillator : public TemporalKernel
{
public:
    DampedOscillator(double phaseDuration, double dampingFactor, double delay);
    ~DampedOscillator();

    // TemporalKernel interface
    double temporal(double t) const;
    complex<double> fourierTransform(double w) const;


private:
    double m_phaseDuration = 1.0;
    double m_dampingFactor = 0.0;
    double m_delay = 0.0;
};

}
lgnSimulator::DampedOscillator createTemporalDampedOscillatorKernel(const YAML::Node &cfg);

#endif // DAMPEDOSCILLATOR_H
