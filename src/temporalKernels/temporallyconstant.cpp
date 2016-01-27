#include "temporallyconstant.h"

using namespace edog;


TemporallyConstant::TemporallyConstant(double constant)
    : m_constant(constant)
{

}

TemporallyConstant::~TemporallyConstant()
{

}



double TemporallyConstant::temporal(double t)
{
    (void)t;
    return m_constant;
}

complex<double> TemporallyConstant::fourierTransform(double w)
{
    return m_constant * Functions::delta(w,0);
}


TemporallyConstant createTemporallyConstantKernel(const YAML::Node *cfg)
{
    double constant = (*cfg)["constant"].as<double>();

    return TemporallyConstant(constant);

}
