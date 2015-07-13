#include "patchgrating.h"

PatchGrating::PatchGrating(const Config *cfg)
    : Stimuli(cfg)
{
    const Setting & root = cfg->getRoot();
    m_contrast = root["stimuliSettings"]["C"];
    m_spotDiameter = root["stimuliSettings"]["d"];
}

PatchGrating::~PatchGrating()
{

}

double PatchGrating::real(vec2 rVec, double t)
{
    double r = sqrt(dot(rVec,rVec));
    double s = (1 - m_math.heaviside(r - m_spotDiameter * 0.5))
            * m_contrast * cos(dot(rVec, m_k) - m_w * t);

    return s;

}

double PatchGrating::complex(vec2 kVec, double w)
{
    double arg = sqrt(dot(kVec - m_k, kVec - m_k));
    double s = m_math.secondKindBesselFunction(arg * m_spotDiameter * 0.5)
            * m_math.delta(w, m_w) * m_spotDiameter * m_contrast/arg;

    return s*2*PI*PI;

}
