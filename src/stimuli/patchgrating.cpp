#include "patchgrating.h"

PatchGrating::PatchGrating(Integrator integrator, vec2 kd,
                           double wd, double contrast,
                           double spotDiameter)
    : Stimuli(integrator)
    , m_k(kd)
    , m_w(wd)
    , m_contrast(contrast)
    , m_spotDiameter(spotDiameter)
{
}

PatchGrating::~PatchGrating()
{

}

double PatchGrating::valueAtPoint(vec2 rVec, double t)
{
    vec dr = rVec;
    double r = sqrt(dot(dr, dr));
    double s = m_contrast * (1 - Functions::heaviside(r - m_spotDiameter * 0.5))
             * cos(dot(m_k, dr) - m_w * t);

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
