#ifndef DECAYINGEXPONENTIAL_H
#define DECAYINGEXPONENTIAL_H

#include "temporalKernels/temporalkernel.h"

class DecayingExponential : public TemporalKernel
{
public:
    DecayingExponential(double tau, double delay);
    ~DecayingExponential();

    // TemporalKernel interface
    double real(double t);
    double complex(double w);


private:
    double m_tau, m_delay = 0.0;
};

#endif // DECAYINGEXPONENTIAL_H
