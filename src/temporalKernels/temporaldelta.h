#ifndef TEMPORALDELTA_H
#define TEMPORALDELTA_H

#include "temporalKernels/temporalkernel.h"


namespace edog {
class TemporalDelta : public TemporalKernel
{
public:
    TemporalDelta(double delay);
    ~TemporalDelta();

    // TemporalKernel interface
public:
    virtual double temporal(double t);
    virtual complex<double> fourierTransform(double w);

private:
    double m_delay = 0;
};

}
edog::TemporalDelta createTemporalDeltaKernel(const YAML::Node *cfg);

#endif // TEMPORALDELTA_H