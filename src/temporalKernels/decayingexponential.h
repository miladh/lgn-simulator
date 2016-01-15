#ifndef DECAYINGEXPONENTIAL_H
#define DECAYINGEXPONENTIAL_H

#include "temporalKernels/temporalkernel.h"

namespace edog {

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

}
edog::DecayingExponential createDecayingExponentialTemporalKernel(const YAML::Node *cfg);

#endif // DECAYINGEXPONENTIAL_H
