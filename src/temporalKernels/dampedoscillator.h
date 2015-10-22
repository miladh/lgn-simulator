#ifndef DAMPEDOSCILLATOR_H
#define DAMPEDOSCILLATOR_H


#include "temporalKernels/temporalkernel.h"

class DampedOscillator : public TemporalKernel
{
public:
    DampedOscillator(double phaseDuration, double weight);
    ~DampedOscillator();

    // TemporalKernel interface
    double temporal(double t);
    double fourierTransform(double w);


private:
    double m_phaseDuration;
    double m_weight;
};

DampedOscillator createDampedOscillatorTemporalKernel(const Config *cfg);

#endif // DAMPEDOSCILLATOR_H
