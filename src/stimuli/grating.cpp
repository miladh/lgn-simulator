#include "grating.h"

Grating::Grating(const Config *cfg, Integrator integrator)
    : Stimuli(cfg, integrator)
{
    const Setting & root = cfg->getRoot();
    m_contrast = root["stimuliSettings"]["C"];
}

Grating::~Grating()
{

}

double Grating::valueAtPoint(vec2 rVec, double t)
{
    vec dr = rVec;
    double s = m_contrast * cos(dot(m_k, dr) - m_w * t);

    return s;
}

double Grating::fourierTransformAtFrequency(vec2 k, double w)
{
    double s = Functions::delta(k[0], m_k[0])
            * Functions::delta(k[1], m_k[1])
            * Functions::delta(w, -m_w);

    return 8*PI*PI*PI* m_contrast * s;
}

