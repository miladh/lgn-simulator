#include "temporaldelta.h"

using namespace lgnSimulator;


TemporalDelta::TemporalDelta(double delay, double temporalResolution)
    : m_delay(delay)
    , m_peak(1./temporalResolution)
{
}

TemporalDelta::~TemporalDelta()
{

}


double TemporalDelta::temporal(double t) const
{
    return m_peak * Special::delta(t, m_delay);
}

complex<double> TemporalDelta::fourierTransform(double w) const
{
    // Note the additional -1 factor that accounts for the difference
    // in definition of temporal and spatial inverse FT.
    return exp(core::i * m_delay * -w);
}


TemporalDelta createTemporalDeltaKernel(const YAML::Node &cfg)
{
    double delay = cfg["delay"].as<double>();
    double dt = cfg["dt"].as<double>();

    return TemporalDelta(delay, dt);


}
