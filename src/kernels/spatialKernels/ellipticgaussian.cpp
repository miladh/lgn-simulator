#include "ellipticgaussian.h"


/*!
  \class lgnSimulator::EllipticGaussian
  \inmodule lgnSimulator
  \ingroup lgnSimulator-spatialKernel
  \brief Spatial elliptic Gaussian kernel.
 */

using namespace lgnSimulator;


EllipticGaussian::EllipticGaussian(double angle, double widthLong, double widthNarrow)
    :m_angle(angle*core::pi/180.)
    , m_widthLong(widthLong)
    , m_widthNarrow(widthNarrow)
    , m_cosTheta(cos(m_angle))
    , m_sinTheta(sin(m_angle))

{

}

EllipticGaussian::~EllipticGaussian()
{

}

double EllipticGaussian::spatial(vec2 r) const
{

    double exp1 = (r[0] * m_cosTheta + r[1] * m_sinTheta)/m_widthLong;
    double exp2 = (r[1] * m_cosTheta - r[0] * m_sinTheta)/m_widthNarrow;
    return 1./core::pi/m_widthNarrow/m_widthLong * exp(-exp1*exp1 - exp2*exp2);
}

complex<double> EllipticGaussian::fourierTransform(vec2 k) const
{
    double exp1 = (k[0] * m_cosTheta + k[1] * m_sinTheta)*m_widthLong/2;
    double exp2 = (k[1] * m_cosTheta - k[0] * m_sinTheta)*m_widthNarrow/2;

    return  exp(-exp1*exp1 - exp2*exp2);
}



EllipticGaussian createSpatialEllipticGaussianKernel(const YAML::Node &cfg)
{
    double angle = cfg["angle"].as<double>();
    double widthLong = cfg["widthLong"].as<double>();
    double widthNarrow = cfg["widthNarrow"].as<double>();

    return EllipticGaussian(angle, widthLong, widthNarrow);
}
