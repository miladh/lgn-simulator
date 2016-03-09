#include "temporallyconstant.h"

using namespace lgnSimulator;


TemporallyConstant::TemporallyConstant(double constant, double temporalFreqResolution)
    : m_constant(constant)
    , m_peak(1./temporalFreqResolution)
{
}

TemporallyConstant::~TemporallyConstant()
{

}


double TemporallyConstant::temporal(double t) const
{
    (void)t;
    return m_constant;
}

complex<double> TemporallyConstant::fourierTransform(double w) const
{
    return m_constant * m_peak * 2. * core::pi * Special::delta(w,0);
}


TemporallyConstant createTemporallyConstantKernel(const YAML::Node &cfg)
{
    double constant = cfg["constant"].as<double>();
    double dt = cfg["dt"].as<double>();
    double nt = cfg["nt"].as<double>();
    double peak = 2.*core::pi/pow(2,nt)/dt;

    return TemporallyConstant(constant, peak);

}
