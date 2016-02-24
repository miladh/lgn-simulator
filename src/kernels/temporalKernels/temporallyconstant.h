#ifndef TEMPORALLYCONSTANT_H
#define TEMPORALLYCONSTANT_H

#include "kernels/temporalKernels/temporalkernel.h"


namespace lgnSimulator {
class TemporallyConstant : public TemporalKernel
{
public:
    TemporallyConstant(double constant, double temporalFreqResolution);
    ~TemporallyConstant();

    // TemporalKernel interface
public:
    double temporal(double t) const;
    virtual complex<double> fourierTransform(double w) const;

private:
    double m_constant= 0.0;
    double m_peak = 1.0;
};

}

lgnSimulator::TemporallyConstant createTemporallyConstantKernel(const YAML::Node &cfg);

#endif // TEMPORALLYCONSTANT_H
