#include "circlemaskgrating.h"

using namespace lgnSimulator;


CircleMaskGrating::CircleMaskGrating(const Integrator &integrator,
                                     double spatialFreq, double orientation,
                                     double temporalFreq,double contrast,
                                     double maskSize)
    : Grating(integrator, spatialFreq, orientation, temporalFreq, contrast, maskSize)
{
    m_mask= "circle";
}

CircleMaskGrating::~CircleMaskGrating()
{

}


double CircleMaskGrating::valueAtPoint(vec2 rVec, double t) const
{
    double r = sqrt(dot(rVec, rVec));

    double s = m_contrast * (1 - Special::heaviside(r  - m_maskSize * 0.5))
            * cos(dot(m_kVec, rVec) - m_w * t);

    return s;
}



complex<double> CircleMaskGrating::fourierTransformAtFrequency(vec2 kVec, double w) const
{
    if(!Special::delta(w, m_w)){
        return 0;
    }

    double s = m_contrast * core::pi * core::pi * m_maskSize * m_maskSize * 0.5;
    double arg = sqrt(dot(kVec - m_kVec, kVec - m_kVec)) * m_maskSize * 0.5;

    if(arg != 0){
        s *= 2. * Special::secondKindBesselFunction(arg)/arg;
    }

    return s/m_integrator.temporalFreqResolution();
}
