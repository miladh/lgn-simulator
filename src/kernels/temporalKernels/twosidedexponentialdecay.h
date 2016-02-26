#ifndef TWOSIDEDEXPONENTIALDECAY_H
#define TWOSIDEDEXPONENTIALDECAY_H

#include "temporalkernel.h"

namespace lgnSimulator {
class TwoSidedExponentialDecay : public TemporalKernel
{
public:
    TwoSidedExponentialDecay(double tau, double delay);

    // TemporalKernel interface
    virtual double temporal(double t) const;
    virtual complex<double> fourierTransform(double w) const;



private:
    double m_tau = 0.0;
    double m_delay = 0.0;
};

}

lgnSimulator::TwoSidedExponentialDecay createTemporalTwoSidedExponentialDecayKernel(const YAML::Node &cfg);
#endif // TWOSIDEDEXPONENTIALDECAY_H
