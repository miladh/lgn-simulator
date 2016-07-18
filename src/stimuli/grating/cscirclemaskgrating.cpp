#include "cscirclemaskgrating.h"


using namespace lgnSimulator;
CSCircleMaskGrating::CSCircleMaskGrating(Integrator* const integrator,
                                         double spatialFreq, double orientation,
                                         double temporalFreq,double contrast,
                                         double surroundSize)
    : Grating(integrator, spatialFreq, orientation, temporalFreq, contrast, surroundSize)
{
   m_mask= "cscircle";
   m_centerSize = 5;
//   m_centerK = vec2{m_k*cos(core::pi/2), m_k*sin(core::pi/2)};
   m_centerK = m_kVec*2;

}

void CSCircleMaskGrating::computeFourierTransform()
{
    computeSpatiotemporal();
    m_fourierTransform = m_integrator->forwardFFT(m_spatiotemporal);

}

CSCircleMaskGrating::~CSCircleMaskGrating()
{

}


double CSCircleMaskGrating::valueAtPoint(vec2 rVec, double t) const
{
    double r = sqrt(dot(rVec, rVec));

    if(r  > m_maskSize * 0.5){
        return 0;
    }else if( r > m_centerSize * 0.5){
        return  m_contrast *  cos(dot(m_kVec, rVec) - m_w * t);
    }else{
        return m_contrast * cos(dot(m_centerK, rVec) - m_w * t);

    }
}



complex<double> CSCircleMaskGrating::fourierTransformAtFrequency(vec2 k, double w) const
{

//    double arg1 = sqrt(dot(k - m_kVec, k - m_kVec))* m_maskSize * 0.5;
//    double arg2 = sqrt(dot(k + m_kVec, k + m_kVec))* m_maskSize * 0.5;


//    double term1 = Special::delta(w, m_w);
//    double term2 = Special::delta(w, -m_w);

//    if(arg1!= 0){
//        term1 *= 2*Special::secondKindBesselFunction(arg1)/arg1;
//    }
//    if(arg2!= 0){
//        term2 *= 2*Special::secondKindBesselFunction(arg2)/arg2;
//    }

//    return m_contrast * core::pi * core::pi * m_maskSize * m_maskSize * 0.25
//            * (term1 + term2) / m_integrator->temporalFreqResolution();

}
