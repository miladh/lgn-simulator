#include "patchgrating.h"

PatchGrating::PatchGrating(const Config *cfg)
    : Stimuli(cfg)
{
    const Setting & root = cfg->getRoot();
    m_C = root["stimuliSettings"]["C"];
    m_d = root["stimuliSettings"]["d"];
}

PatchGrating::~PatchGrating()
{

}

double PatchGrating::real(vec2 rVec, double t)
{
    double r = sqrt(dot(rVec,rVec));
    double s = m_C*(1 - m_math.heaviside(r - m_d * 0.5))
            * cos(dot(rVec, m_k) - m_w * t);

    return s;
}

double PatchGrating::complex(vec2 kVec, double w)
{
    double arg = sqrt(dot(kVec - m_k, kVec - m_k));
    double s = m_math.secondKindBesselFunction(arg * m_d * 0.5)
            * m_math.delta(w, m_w);
    return s*2*PI*PI*m_d*m_C/arg;

}
