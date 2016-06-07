#include "dog.h"


using namespace lgnSimulator;


DOG::DOG(double a, double b, double c):
    m_relativeStrength(c)
{
    m_centre = new SpatialGaussian(a);
    m_surround = new SpatialGaussian(b);
}

DOG::~DOG()
{
}

double DOG::spatial(vec2 r) const
{
    return m_centre->spatial(r) - m_relativeStrength * m_surround->spatial(r);
}

complex<double> DOG::fourierTransform(vec2 k) const
{
    return m_centre->fourierTransform(k)
            - m_relativeStrength * m_surround->fourierTransform(k);
}



DOG createSpatialDOGKernel(const YAML::Node &cfg)
{

    double a = cfg["a"].as<double>();
    double b = cfg["b"].as<double>();
    double c = cfg["c"].as<double>();

    return DOG(a, b, c);
}
