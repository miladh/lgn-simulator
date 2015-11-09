#ifndef TEMPORALDELTA_H
#define TEMPORALDELTA_H

#include "temporalkernel.h"



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

TemporalDelta createTemporalDeltaKernel(const Config *cfg);

#endif // TEMPORALDELTA_H
