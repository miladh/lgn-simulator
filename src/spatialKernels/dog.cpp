#include "dog.h"


using namespace lgnSimulator;


DOG::DOG(double A, double a, double B, double b)
{
    m_centre = new Gaussian(A,a);
    m_surround = new Gaussian(B,b);
}

DOG::~DOG()
{
}

double DOG::spatial(vec2 r)
{
    return m_centre->spatial(r) - m_surround->spatial(r);
}

complex<double> DOG::fourierTransform(vec2 k)
{
    return m_centre->fourierTransform(k) - m_surround->fourierTransform(k);
}



DOG createDOGSpatialKernel(const YAML::Node *cfg)
{

    double A = (*cfg)["A"].as<double>();
    double a = (*cfg)["a"].as<double>();
    double B = (*cfg)["B"].as<double>();
    double b = (*cfg)["b"].as<double>();


    return DOG(A, a, B, b);
}
