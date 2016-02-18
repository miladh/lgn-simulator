#ifndef DAMPEDOSCILLATOR_H
#define DAMPEDOSCILLATOR_H


#include "kernels/temporalKernels/temporalkernel.h"

namespace lgnSimulator {

class DampedOscillator : public TemporalKernel
{
public:
    DampedOscillator(double phaseDuration, double weight);
    ~DampedOscillator();

    // TemporalKernel interface
    double temporal(double t) const;
    complex<double> fourierTransform(double w) const;


private:
    double m_phaseDuration;
    double m_weight;
};

}
lgnSimulator::DampedOscillator createTemporalDampedOscillatorKernel(const YAML::Node &cfg);

#endif // DAMPEDOSCILLATOR_H