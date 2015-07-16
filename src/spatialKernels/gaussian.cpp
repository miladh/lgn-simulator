#include "gaussian.h"

Gaussian::Gaussian(double weight, double spread)
    : m_weight(weight)
    , m_spread(spread)
{

}

Gaussian::~Gaussian()
{

}

double Gaussian::real(vec2 rVec)
{
    double r = sqrt(dot(rVec,rVec));
    throw "not implemented";
}

double Gaussian::complex(vec2 kVec)
{
    double k = sqrt(dot(kVec,kVec));
    return m_weight * exp(-k*k * m_spread*m_spread * 0.25);
}

