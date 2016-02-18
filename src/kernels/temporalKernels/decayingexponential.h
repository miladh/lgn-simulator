#ifndef DECAYINGEXPONENTIAL_H
#define DECAYINGEXPONENTIAL_H

#include "kernels/temporalKernels/temporalkernel.h"

namespace lgnSimulator {

class DecayingExponential : public TemporalKernel
{
public:
    DecayingExponential(double tau, double delay);
    ~DecayingExponential();

    // TemporalKernel interface
    double temporal(double t) const;
    complex<double> fourierTransform(double w) const;


private:
    double m_tau = 0.0;
    double m_delay = 0.0;
};

}
lgnSimulator::DecayingExponential createTemporalDecayingExponentialKernel(const YAML::Node &cfg);

#endif // DECAYINGEXPONENTIAL_H
