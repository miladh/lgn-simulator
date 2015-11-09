#include "temporaldelta.h"

TemporalDelta::TemporalDelta(double delay)
    :m_delay(delay)
{

}

TemporalDelta::~TemporalDelta()
{

}


double TemporalDelta::temporal(double t)
{
    return Functions::delta(t, m_delay);
}

complex<double> TemporalDelta::fourierTransform(double w)
{
    return cos(m_delay * w);
}


TemporalDelta createTemporalDeltaKernel(const Config *cfg)
{
    const Setting & root = cfg->getRoot();
    double delay = root["temporalKernelSettings"]["delay"];

    return TemporalDelta(delay);
}
