#include "temporallyconstant.h"

/*!
  \class lgnSimulator::TemporallyConstant
  \inmodule lgnSimulator
  \ingroup lgnSimulator-temporalKernel
  \brief Constant temporal kernel.
 */

using namespace lgnSimulator;


TemporallyConstant::TemporallyConstant(double temporalInterval,
                                       double temporalFreqResolution)
    : m_temporalInterval(temporalInterval)
    , m_peak(1./temporalFreqResolution)
{
}

TemporallyConstant::~TemporallyConstant()
{

}


double TemporallyConstant::temporal(double t) const
{
    (void)t;
    return 1./m_temporalInterval;
}

complex<double> TemporallyConstant::fourierTransform(double w) const
{
    return 1./m_temporalInterval * m_peak * 2. * core::pi * Special::delta(w,0);
}


TemporallyConstant createTemporallyConstantKernel(const YAML::Node &cfg)
{
    double dt = cfg["dt"].as<double>();
    double nt = cfg["nt"].as<double>();

    double temporalInterval = pow(2. ,nt)*dt;
    double peak = 2.*core::pi/temporalInterval;

    return TemporallyConstant(temporalInterval, peak);

}
