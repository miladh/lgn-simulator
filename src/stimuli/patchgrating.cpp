#include "patchgrating.h"

PatchGrating::PatchGrating(const Config *cfg, Integrator integrator)
    : Stimuli(cfg, integrator)
{
    const Setting & root = cfg->getRoot();
    m_contrast = root["stimuliSettings"]["C"];
    m_spotDiameter = root["stimuliSettings"]["d"];
}

PatchGrating::~PatchGrating()
{

}

double PatchGrating::valueAtPoint(vec2 rVec, double t)
{
    vec dr = rVec - vec{0.5, 0.5};
    double r = sqrt(dot(dr, dr));
    double s = (1 - Functions::heaviside(r - m_spotDiameter * 0.5))
            * m_contrast * cos(dot(m_k, dr) - m_w * t);

    return s;

}

double PatchGrating::fourierTransformAtFrequency(vec2 kVec, double w)
{
    double arg = sqrt(dot(kVec - m_k, kVec - m_k)) * m_spotDiameter * 0.5;
    double s = m_contrast * PI * PI * m_spotDiameter * m_spotDiameter * 0.5
            * Functions::delta(w, -m_w);
    if(arg != 0){
        s *= 2. * Functions::secondKindBesselFunction(arg)/arg;
    }
    return s;

}


//double arg1 = sqrt(dot(kVec - m_k, kVec - m_k)) * m_spotDiameter * 0.5;
//double arg2 = sqrt(dot(kVec + m_k, kVec + m_k)) * m_spotDiameter * 0.5;
//double s1 = Functions::delta(w, m_w);
//double s2 = Functions::delta(w, -m_w);

//if(arg1 != 0 && s1!=0){
//    s1 *= 2. * Functions::secondKindBesselFunction(arg1)/arg1;
//}

//if(arg2 != 0 && s2!=0){
//    s2 *= 2. * Functions::secondKindBesselFunction(arg2)/arg2;
//}
//return (s1+s2)*0.5*m_contrast * PI * PI * m_spotDiameter * m_spotDiameter * 0.5;
