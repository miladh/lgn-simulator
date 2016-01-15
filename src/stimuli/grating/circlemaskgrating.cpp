#include "circlemaskgrating.h"

using namespace edog;


CircleMaskGrating::CircleMaskGrating(Integrator *integrator,
                           vec2 kd, double wd, double contrast, double maskSize)
    : Grating(integrator, kd, wd, contrast, maskSize)
{
}

CircleMaskGrating::~CircleMaskGrating()
{

}

double CircleMaskGrating::valueAtPoint(vec2 rVec, double t)
{
    double r = sqrt(dot(rVec, rVec));
    double s = m_contrast * (1 - Functions::heaviside(r - m_maskSize * 0.5))
            * cos(dot(m_k, rVec) - m_w * t);

    return s;
}



double CircleMaskGrating::fourierTransformAtFrequency(vec2 kVec, double w)
{
    if(!Functions::delta(w, -m_w)){
        return 0;
    }

    double s = m_contrast * PI * PI * m_maskSize * m_maskSize * 0.5;
    double arg = sqrt(dot(kVec - m_k, kVec - m_k)) * m_maskSize * 0.5;

    if(arg != 0){
        s *= 2. * Functions::secondKindBesselFunction(arg)/arg;
    }

    return s/m_integrator->temporalFreqResolution();
}
