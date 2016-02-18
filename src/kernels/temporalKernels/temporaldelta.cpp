#include "temporaldelta.h"

using namespace lgnSimulator;


TemporalDelta::TemporalDelta(double delay)
    :m_delay(delay)
{

}

TemporalDelta::~TemporalDelta()
{

}


double TemporalDelta::temporal(double t) const
{
    return SpecialFunctions::delta(t, m_delay);
}

complex<double> TemporalDelta::fourierTransform(double w) const
{
    // Note that we have +i, and not -i due to the definition of
    // temporal FT.
    return exp(core::i * m_delay * w);
}


TemporalDelta createTemporalDeltaKernel(const YAML::Node &cfg)
{
    double delay = cfg["delay"].as<double>();

    return TemporalDelta(delay);
}
