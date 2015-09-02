#include "dogstim.h"

DOGstim::DOGstim(const Config *cfg)
    : Stimuli(cfg)
{

}

DOGstim::~DOGstim()
{

}

double DOGstim::spatial(vec2 rVec, double t)
{
    return m_dog.real(rVec);

}

double DOGstim::frequency(vec2 k, double w)
{
    m_dog.complex(k);
}

