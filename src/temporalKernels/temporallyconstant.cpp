#include "temporallyconstant.h"

TemporallyConstant::TemporallyConstant(double constant)
    : m_constant(constant)
{

}

TemporallyConstant::~TemporallyConstant()
{

}



double TemporallyConstant::temporal(double t)
{
    return m_constant;
}

double TemporallyConstant::fourierTransform(double w)
{
    return m_constant * Functions::delta(w,0);
}


TemporallyConstant createTemporallyConstantTemporalKernel(const Config *cfg)
{
    const Setting & root = cfg->getRoot();
    double constant = root["temporalKernelSettings"]["constant"];

    return TemporallyConstant(constant);

}
