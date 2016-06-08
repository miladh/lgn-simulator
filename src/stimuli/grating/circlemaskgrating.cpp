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

//    cout << (1 - Special::heaviside(r  - m_maskSize * 0.5)) << endl;
    double s = m_contrast * (1 - Special::heaviside(r  - m_maskSize * 0.5))
            * cos(dot(m_kVec, rVec) - m_w * t);

    return s;
}


complex<double> CircleMaskGrating::fourierTransformAtFrequency(vec2 k, double w) const
{

    double arg1 = sqrt(dot(k - m_kVec, k - m_kVec))* m_maskSize * 0.5;
    double arg2 = sqrt(dot(k + m_kVec, k + m_kVec))* m_maskSize * 0.5;


    double term1 = Special::delta(w, m_w);
    double term2 = Special::delta(w, -m_w);

    if(Special::delta(k, m_kVec) != 1){
        term1*=2.*Special::secondKindBesselFunction(arg1)/arg1;
    }
    if(Special::delta(k, -m_kVec)!= 1){
        term2*=2.*Special::secondKindBesselFunction(arg2)/arg2;
    }

    return m_contrast * core::pi * core::pi * m_maskSize * m_maskSize * 0.5
           * 0.5 * (term1 + term2)
            /m_integrator.temporalFreqResolution();
//            / m_integrator.spatialFreqResolution()
//            / m_integrator.spatialFreqResolution();
}


//    double s = (Special::delta(k, m_kVec) * Special::delta(w, m_w)
//                + Special::delta(k, -m_kVec) * Special::delta(w, -m_w));

//    return 4.*core::pi*core::pi*core::pi * m_contrast
//            / m_integrator.temporalFreqResolution()
//            / m_integrator.spatialFreqResolution()
//            / m_integrator.spatialFreqResolution() * s;

//double s = m_contrast * core::pi * core::pi * m_maskSize * m_maskSize * 0.5;
//double arg = sqrt(dot(k - m_kVec, k - m_kVec)) * m_maskSize * 0.5;

//if(arg != 0){
//    s *= 2. * Special::secondKindBesselFunction(arg)/arg;
//}

//retur/*n s/m_integrator.temporalFreqResolution();
//}*/
