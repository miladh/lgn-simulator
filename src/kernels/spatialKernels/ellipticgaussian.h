#ifndef ELLIPTICGAUSSIAN_H
#define ELLIPTICGAUSSIAN_H

#include "spatialkernel.h"
namespace lgnSimulator {
class EllipticGaussian : public SpatialKernel
{
public:
    EllipticGaussian(double angle,double widthLong, double widthNarrow);
    ~EllipticGaussian();

    // SpatialKernel interface
    double spatial(vec2 r) const;
    complex<double> fourierTransform(vec2 k) const;

private:
    double m_angle = 0.0;
    double m_widthLong = 0.0;
    double m_widthNarrow = 0.0;
    double m_cosTheta = 0.0;
    double m_sinTheta = 0.0;
};

}

lgnSimulator::EllipticGaussian createSpatialEllipticGaussianKernel(const YAML::Node &cfg);

#endif // ELLIPTICGAUSSIAN_H
