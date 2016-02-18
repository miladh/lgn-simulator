#ifndef TEMPORALDELTA_H
#define TEMPORALDELTA_H

#include "kernels/temporalKernels/temporalkernel.h"


namespace lgnSimulator {
class TemporalDelta : public TemporalKernel
{
public:
    TemporalDelta(double delay);
    ~TemporalDelta();

    // TemporalKernel interface
public:
    virtual double temporal(double t) const;
    virtual complex<double> fourierTransform(double w) const;

private:
    double m_delay = 0;
};

}
lgnSimulator::TemporalDelta createTemporalDeltaKernel(const YAML::Node &cfg);

#endif // TEMPORALDELTA_H
