#ifndef DOE_H
#define DOE_H

#include "temporalkernel.h"

namespace lgnSimulator {

class DOE : public TemporalKernel
{
public:
    DOE(double cenLatency, double surLatency, double delay);

    // TemporalKernel interface
    virtual double temporal(double t) const;
    virtual complex<double> fourierTransform(double w) const;

private:
    double m_cenLatency = 0.0;
    double m_surLatency = 0.0;
    double m_delay = 0.0;
};

}

lgnSimulator::DOE createTemporalDOEKernel(const YAML::Node &cfg);
#endif // DOE_H
