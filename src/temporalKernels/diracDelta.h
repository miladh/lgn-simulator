#ifndef LINEAR_H
#define LINEAR_H

#include "temporalkernel.h"
#include "math/functions.h"

class DiracDelta : public TemporalKernel
{
public:
    DiracDelta(double t);
    ~DiracDelta();

    // TemporalKernel interface
public:
    double temporal(double t);
    double fourierTransform(double w);

private:
    double m_t = 0.;
};

#endif // LINEAR_H
