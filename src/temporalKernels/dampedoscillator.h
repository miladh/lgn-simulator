#ifndef DAMPEDOSCILLATOR_H
#define DAMPEDOSCILLATOR_H


#include "temporalKernels/temporalkernel.h"

namespace lgnSimulator {

class DampedOscillator : public TemporalKernel
{
public:
    DampedOscillator(double phaseDuration, double weight);
    ~DampedOscillator();

    // TemporalKernel interface
    double temporal(double t);
    complex<double> fourierTransform(double w);


private:
    double m_phaseDuration;
    double m_weight;
};

}
lgnSimulator::DampedOscillator createTemporalDampedOscillatorKernel(const YAML::Node *cfg);

#endif // DAMPEDOSCILLATOR_H
