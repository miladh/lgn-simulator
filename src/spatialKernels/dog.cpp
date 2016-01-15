#include "dog.h"

using namespace edog;


DOG::DOG(double A, double a, double B, double b)
    : m_A(A)
    , m_a(a)
    , m_B(B)
    , m_b(b)
{

}

DOG::~DOG()
{
}

double DOG::spatial(vec2 r)
{

    double r2 = dot(r,r);
    double center   = m_A / (m_a*m_a) / PI * exp(-r2 / (m_a*m_a));
    double surround = m_B / (m_b*m_b) / PI * exp(-r2 / (m_b*m_b));


    return center - surround;
}

double DOG::fourierTransform(vec2 k)
{
    double k2 = dot(k,k);
    double center   = m_A * exp(-k2 * m_a*m_a / 4.);
    double surround = m_B * exp(-k2 * m_b*m_b / 4.);


    return center - surround;
}



DOG createDOGSpatialKernel(const YAML::Node *cfg)
{

    double A = (*cfg)["A"].as<double>();
    double a = (*cfg)["a"].as<double>();
    double B = (*cfg)["B"].as<double>();
    double b = (*cfg)["b"].as<double>();


    return DOG(A, a, B, b);
}
