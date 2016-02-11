#ifndef TEMPORALGAUSSIAN_H
#define TEMPORALGAUSSIAN_H

#include "temporalkernel.h"


namespace lgnSimulator {
class TemporalGaussian : public TemporalKernel
{
public:
    TemporalGaussian(double A, double a, double delay);

    virtual double temporal(double t);
    virtual complex<double> fourierTransform(double w);


private:
    double m_A = 0.0;
    double m_a = 0.0;
    double m_delay = 0.0;
};

}

lgnSimulator::TemporalGaussian createGaussianTemporalKernel(const YAML::Node *cfg);

#endif // TEMPORALGAUSSIAN_H
