#ifndef DAMPEDOSCILLATOR_H
#define DAMPEDOSCILLATOR_H


#include "temporalKernels/temporalkernel.h"

class DampedOscillator : public TemporalKernel
{
public:
    DampedOscillator(double phaseDuration, double weight);
    ~DampedOscillator();

    // TemporalKernel interface
    double real(double t);
    double complex(double w);


private:
    double m_phaseDuration;
    double m_weight;
};

#endif // DAMPEDOSCILLATOR_H
