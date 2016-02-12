#include "ellipticgaussian.h"

using namespace lgnSimulator;


EllipticGaussian::EllipticGaussian(double weight, double angle,
                                   double widthLong, double widthNarrow)
    : m_weight(weight)
    , m_angle(angle*PI/180.)
    , m_widthLong(widthLong)
    , m_widthNarrow(widthNarrow)
    , m_cosTheta(cos(m_angle))
    , m_sinTheta(sin(m_angle))

{

}

EllipticGaussian::~EllipticGaussian()
{

}

double EllipticGaussian::spatial(vec2 r)
{

    double exp1 = (r[0] * m_cosTheta + r[1] * m_sinTheta)/m_widthLong;
    double exp2 = (r[1] * m_cosTheta - r[0] * m_sinTheta)/m_widthNarrow;
    return m_weight/PI/m_widthNarrow/m_widthLong * exp(-exp1*exp1 - exp2*exp2);
}

complex<double> EllipticGaussian::fourierTransform(vec2 k)
{
    double exp1 = (k[0] * m_cosTheta + k[1] * m_sinTheta)*m_widthLong/2;
    double exp2 = (k[1] * m_cosTheta - k[0] * m_sinTheta)*m_widthNarrow/2;

    return m_weight * exp(-exp1*exp1 - exp2*exp2);
}



EllipticGaussian createSpatialEllipticGaussianKernel(const YAML::Node *cfg)
{
    double weight = (*cfg)["weight_elliptic"].as<double>();
    double angle = (*cfg)["angle"].as<double>();
    double widthLong = (*cfg)["widthLong"].as<double>();
    double widthNarrow = (*cfg)["widthNarrow"].as<double>();

    return EllipticGaussian(weight, angle, widthLong, widthNarrow);
}
