#include "dog.h"

/*!
  \class lgnSimulator::DOG
  \inmodule lgnSimulator
  \ingroup lgnSimulator-spatialKernel
  \brief Spatial difference of Gaussians kernel.
 */


using namespace lgnSimulator;


DOG::DOG(double a, double b, double c)
    : m_relativeStrength(c)
    , m_center(new SpatialGaussian(a))
    , m_surround( new SpatialGaussian(b))
{
}

DOG::DOG(const DOG &dog)
{
    m_center = new SpatialGaussian(1);
    *m_center = *dog.m_center;

    m_surround = new SpatialGaussian(1);
    *m_surround = *dog.m_surround;

    m_relativeStrength = dog.m_relativeStrength;
}

DOG::~DOG()
{
    delete m_center;
    delete m_surround;
}

double DOG::spatial(vec2 r) const
{
    return m_center->spatial(r) - m_relativeStrength * m_surround->spatial(r);
}

complex<double> DOG::fourierTransform(vec2 k) const
{
    return m_center->fourierTransform(k)
            - m_relativeStrength * m_surround->fourierTransform(k);
}

double DOG::relativeStrength() const
{
    return m_relativeStrength;
}



DOG createSpatialDOGKernel(const YAML::Node &cfg)
{

    double a = cfg["a"].as<double>();
    double b = cfg["b"].as<double>();
    double c = cfg["c"].as<double>();

    return DOG(a, b, c);
}
