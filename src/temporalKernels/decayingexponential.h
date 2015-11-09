#ifndef DECAYINGEXPONENTIAL_H
#define DECAYINGEXPONENTIAL_H

#include "temporalKernels/temporalkernel.h"

class DecayingExponential : public TemporalKernel
{
public:
    DecayingExponential(double tau, double delay);
    ~DecayingExponential();

    // TemporalKernel interface
    double temporal(double t);
    complex<double> fourierTransform(double w);


private:
    double m_tau = 0.0;
    double m_delay = 0.0;
};

DecayingExponential createDecayingExponentialTemporalKernel(const Config *cfg);

#endif // DECAYINGEXPONENTIAL_H
