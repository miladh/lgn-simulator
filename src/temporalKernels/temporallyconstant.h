#ifndef TEMPORALLYCONSTANT_H
#define TEMPORALLYCONSTANT_H

#include "temporalKernels/temporalkernel.h"


namespace edog {
class TemporallyConstant : public TemporalKernel
{
public:
    TemporallyConstant(double constant);
    ~TemporallyConstant();

    // TemporalKernel interface
public:
    double temporal(double t);
    virtual complex<double> fourierTransform(double w);

private:
    double m_constant= 0.0;
};

}

edog::TemporallyConstant createTemporallyConstantTemporalKernel(const YAML::Node *cfg);

#endif // TEMPORALLYCONSTANT_H
