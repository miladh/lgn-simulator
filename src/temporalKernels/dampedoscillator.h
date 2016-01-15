#ifndef DAMPEDOSCILLATOR_H
#define DAMPEDOSCILLATOR_H


#include "temporalKernels/temporalkernel.h"

namespace edog {

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
edog::DampedOscillator createDampedOscillatorTemporalKernel(const YAML::Node *cfg);

#endif // DAMPEDOSCILLATOR_H
