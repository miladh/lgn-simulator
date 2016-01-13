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

complex<double> TemporallyConstant::fourierTransform(double w)
{
    return m_constant * Functions::delta(w,0);
}


TemporallyConstant createTemporallyConstantTemporalKernel(const YAML::Node *cfg)
{
    double constant = (*cfg)["temporalKernelSettings"]["constant"].as<double>();

    return TemporallyConstant(constant);

}
