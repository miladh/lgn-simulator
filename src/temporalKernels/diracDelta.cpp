#include "diracDelta.h"

DiracDelta::DiracDelta(double t)
    : m_t(t)
{

}

DiracDelta::~DiracDelta()
{

}

double DiracDelta::temporal(double t)
{

    return Functions::delta(t,m_t);
}

double DiracDelta::fourierTransform(double w)
{
    return cos(w * m_t);
}

