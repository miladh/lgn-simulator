#ifndef SPATIALDELTA_H
#define SPATIALDELTA_H

#include "spatialkernel.h"


namespace edog {

class SpatialDelta : public SpatialKernel
{
public:
    SpatialDelta(double weight, vec2 r0);

    // SpatialKernel interface
    virtual double spatial(vec2 rVec);
    virtual complex<double> fourierTransform(vec2 kVec);


private:
    double m_weight = 0.0;
    vec2 m_r0;
};

}

edog::SpatialDelta createSpatialDeltaKernel(const YAML::Node *cfg);

#endif // SPATIALDELTA_H
