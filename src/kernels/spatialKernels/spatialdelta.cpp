#include "spatialdelta.h"

/*!
  \class lgnSimulator::SpatialDelta
  \inmodule lgnSimulator
  \ingroup lgnSimulator-spatialKernel
  \brief Spatial delta kernel.
 */


using namespace lgnSimulator;
SpatialDelta::SpatialDelta(double spatialResolution, vec2 shift)
    : m_peak(1./spatialResolution/spatialResolution)
    , m_shift(shift)
{
}


double lgnSimulator::SpatialDelta::spatial(vec2 r) const
{
    return m_peak * Special::delta(m_shift, r);
}

complex<double> lgnSimulator::SpatialDelta::fourierTransform(vec2 k) const
{
    return exp(-core::i * dot(k, m_shift));
}



SpatialDelta createSpatialDeltaKernel(const YAML::Node &cfg)
{

    double ds = cfg["ds"].as<double>();
    vec2 shift = {cfg["shift"][0].as<double>(),
                  cfg["shift"][1].as<double>()};


    return SpatialDelta(ds, shift);
}
