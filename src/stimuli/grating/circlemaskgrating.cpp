#include "circlemaskgrating.h"

using namespace lgnSimulator;


CircleMaskGrating::CircleMaskGrating(Integrator* const integrator,
                                     double spatialFreq, double temporalFreq,
                                     double contrast, double phase,
                                     double orientation, double maskSize)
    : Grating(integrator, spatialFreq, temporalFreq, contrast, phase, orientation)
    , m_maskSize(maskSize)
{
    m_mask= "circle";

    if(m_maskSize > m_spatialVec.max()-m_spatialVec.min()){
        cerr << "Warning: mask size (" << maskSize
             << ") larger than grid length: " << m_spatialVec.max()-m_spatialVec.min()
             << endl;
    }
}

CircleMaskGrating::~CircleMaskGrating()
{
}


double CircleMaskGrating::valueAtPoint(vec2 rVec, double t) const
{
    double r = sqrt(dot(rVec, rVec));
    double s = m_contrast * (1 - Special::heaviside(r  - m_maskSize * 0.5))
            * cos(dot(m_kVec, rVec) - m_w * t+ m_phase);

    return s;
}

complex<double> CircleMaskGrating::fourierTransformAtFrequency(vec2 k, double w) const
{

    double arg1 = sqrt(dot(k - m_kVec, k - m_kVec))* m_maskSize * 0.5;
    double arg2 = sqrt(dot(k + m_kVec, k + m_kVec))* m_maskSize * 0.5;

    complex<double> term1 = Special::delta(w, m_w) * m_phase_p_ft; //CHECK SIGN!!!
    complex<double> term2 = Special::delta(w, -m_w) * m_phase_m_ft;

    if(arg1!= 0){
        term1 *= 2*Special::secondKindBessel(arg1)/arg1;
    }
    if(arg2!= 0){
        term2 *= 2*Special::secondKindBessel(arg2)/arg2;
    }

    return m_contrast * core::pi * core::pi * m_maskSize * m_maskSize * 0.25
            * (term1 + term2) / m_integrator->temporalFreqResolution();

}

double CircleMaskGrating::maskSize() const
{
    return m_maskSize;
}
