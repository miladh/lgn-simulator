#include "circlemaskgrating.h"

using namespace lgnSimulator;


CircleMaskGrating::CircleMaskGrating(Integrator* const integrator,
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

//void CircleMaskGrating::computeFourierTransform()
//{
//    computeSpatiotemporal();
//    m_fourierTransform = m_integrator->forwardFFT(m_spatiotemporal);

//}



double CircleMaskGrating::valueAtPoint(vec2 rVec, double t) const
{
    double r = sqrt(dot(rVec, rVec));
    double s = m_contrast * (1 - Special::heaviside(r  - m_maskSize * 0.5))
            * cos(dot(m_kVec, rVec) - m_w * t);

    return s;
}

complex<double> CircleMaskGrating::fourierTransformAtFrequency(vec2 k, double w) const
{

    double arg1 = sqrt(dot(k - m_kVec, k - m_kVec))* m_maskSize * 0.5;
    double arg2 = sqrt(dot(k + m_kVec, k + m_kVec))* m_maskSize * 0.5;


    double term1 = Special::delta(w, m_w); // CHECK SIGN!!!
    double term2 = Special::delta(w, -m_w);

    if(arg1!= 0){
        term1 *= 2*Special::secondKindBessel(arg1)/arg1;
    }
    if(arg2!= 0){
        term2 *= 2*Special::secondKindBessel(arg2)/arg2;
    }

    return m_contrast * core::pi * core::pi * m_maskSize * m_maskSize * 0.25
            * (term1 + term2) / m_integrator->temporalFreqResolution();

}
