#include "dogstim.h"

DOGstim::DOGstim(const Config *cfg, Integrator integrator)
    : Stimuli(cfg,integrator)
{

}

DOGstim::~DOGstim()
{

}

double DOGstim::valueAtPoint(vec2 rVec, double t)
{
    vec dr = rVec - vec{0.5, 0.5};
    return m_dog.real(dr) * cos(-m_w * t);

}

double DOGstim::fourierTransformAtFrequency(vec2 k, double w)
{
    return m_dog.complex(k) * Functions::delta(w, -m_w) * 2 * PI;
}

