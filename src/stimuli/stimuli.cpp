#include "stimuli.h"


Stimuli::Stimuli(const Config *cfg)
{
    const Setting & root = cfg->getRoot();
    m_w = root["stimuliSettings"]["w"];
    m_k[0] = root["stimuliSettings"]["kx"];
    m_k[1] = root["stimuliSettings"]["ky"];

}

Stimuli::~Stimuli()
{

}

double Stimuli::w() const
{
    return m_w;
}


mat Stimuli::real() const
{
    return m_real;
}

void Stimuli::setReal(const mat &real)
{
    m_real = real;
}
mat Stimuli::complex() const
{
    return m_complex;
}

void Stimuli::setComplex(const mat &complex)
{
    m_complex = complex;
}
