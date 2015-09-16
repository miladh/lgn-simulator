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

double PatchGrating::spatial(vec2 rVec, double t)
{
    double r = sqrt(dot(rVec,rVec));
    double s = (1 - Functions::heaviside(r - m_spotDiameter * 0.5))
            * m_contrast * cos(dot(rVec, m_k) - m_w * t);

    return s;

}

double PatchGrating::frequency(vec2 kVec, double w)
{
    double arg = sqrt(dot(kVec - m_k, kVec - m_k)) * m_spotDiameter * 0.5;
    double s = m_contrast * PI * PI * m_spotDiameter * m_spotDiameter * 0.5
            /** Functions::delta(w, m_w)*/;
    if(arg != 0){
        s *= 2. * Functions::secondKindBesselFunction(arg)/arg;
    }

    return s;

}
