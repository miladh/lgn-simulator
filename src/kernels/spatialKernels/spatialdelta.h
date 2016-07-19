#ifndef SPATIALDELTA_H
#define SPATIALDELTA_H

#include "spatialkernel.h"


namespace lgnSimulator {

class SpatialDelta : public SpatialKernel
{
public:
    SpatialDelta(double spatialResolution, vec2 shift);

    // SpatialKernel interface
    virtual double spatial(vec2 r) const;
    virtual complex<double> fourierTransform(vec2 k) const;


private:
    double m_peak = 1.0;
    vec2 m_shift = {0,0};
};

}

lgnSimulator::SpatialDelta createSpatialDeltaKernel(const YAML::Node &cfg);

#endif // SPATIALDELTA_H
