#ifndef TEMPORALLYCONSTANT_H
#define TEMPORALLYCONSTANT_H

#include "temporalkernel.h"



class TemporallyConstant : public TemporalKernel
{
public:
    TemporallyConstant(double constant);
    ~TemporallyConstant();

    // TemporalKernel interface
public:
    virtual double temporal(double t);
    virtual double fourierTransform(double w);

private:
    double m_constant= 0.0;
};

TemporallyConstant createTemporallyConstantTemporalKernel(const Config *cfg);

#endif // TEMPORALLYCONSTANT_H
