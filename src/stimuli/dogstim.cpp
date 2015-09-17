#include "dogstim.h"

DOGstim::DOGstim(const Config *cfg)
    : Stimuli(cfg)
{

}

DOGstim::~DOGstim()
{

}

double DOGstim::valueAtPoint(vec2 rVec, double t)
{
    return m_dog.real(rVec);

}

double DOGstim::fourierTransformAtFrequency(vec2 k, double w)
{
    m_dog.complex(k);
}

