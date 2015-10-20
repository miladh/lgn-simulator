#ifndef ELLIPTICGAUSSIAN_H
#define ELLIPTICGAUSSIAN_H

#include "spatialkernel.h"

class EllipticGaussian : public SpatialKernel
{
public:
    EllipticGaussian(double weight, double angle,
                     double widthLong, double widthNarrow);
    ~EllipticGaussian();

    // SpatialKernel interface
    double spatial(vec2 rVec);
    double fourierTransform(vec2 kVec);

private:
    double m_weight = 0.0;
    double m_angle = 0.0;
    double m_widthLong = 0.0;
    double m_widthNarrow = 0.0;
    double m_cosTheta = 0.0;
    double m_sinTheta = 0.0;
};

#endif // ELLIPTICGAUSSIAN_H
