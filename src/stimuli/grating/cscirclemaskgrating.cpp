#include "cscirclemaskgrating.h"


using namespace lgnSimulator;
CSCircleMaskGrating::CSCircleMaskGrating(Integrator* const integrator,
                                         double spatialFreq, double temporalFreq,
                                         double contrast, double phase,
                                         double orientation, double maskSize,
                                         double surroundSpatialFreq, double surroundTemporalFreq,
                                         double surroundContrast, double surroundPhase,
                                         double surroundOrientation, double surroundMaskSize)
    : Grating(integrator, spatialFreq, temporalFreq, contrast, phase, orientation)
    , m_surroundk(surroundSpatialFreq)
    , m_surroundw(surroundTemporalFreq)
    , m_surroundContrast(surroundContrast)
    , m_surroundPhase(surroundPhase*core::pi/180.)
    , m_surroundOrientation(surroundOrientation*core::pi/180.)
    , m_surroundMaskSize(surroundMaskSize)
    , m_maskSize(maskSize)
{
    m_mask= "cscircle";
    m_surroundkVec = {Special::nearestValue(m_spatialFreqs, m_surroundk*cos(m_surroundOrientation)),
                      Special::nearestValue(m_spatialFreqs,  m_surroundk*sin(m_surroundOrientation))};

    setSurroundSpatialFreq(sqrt(dot(m_kVec, m_kVec)));
    setSurroundOrientation(atan2(m_kVec(1), m_kVec(0)));


    if(m_maskSize > m_spatialVec.max()-m_spatialVec.min()){
        cerr << "Warning: mask size (" << maskSize
             << ") larger than grid length: " << m_spatialVec.max()-m_spatialVec.min()
             << endl;
    }

    if(m_surroundMaskSize > m_spatialVec.max()-m_spatialVec.min()){
        cerr << "Warning: surround mask size (" << m_surroundMaskSize
             << ") larger than grid length: " << m_spatialVec.max()-m_spatialVec.min()
             << endl;
    }

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

    if(r  > m_surroundMaskSize * 0.5){
        return 0;
    }else if( r > m_maskSize * 0.5){
        return  m_surroundContrast *cos(dot(m_surroundkVec, rVec) - m_surroundw*t + m_surroundPhase);
    }else{
        return m_contrast * cos(dot(m_kVec, rVec) - m_w * t + m_phase);

    }
}


complex<double> CSCircleMaskGrating::fourierTransformAtFrequency(vec2 k, double w) const
{

    return centerFourierTransformAtFrequency(k,w)+surroundFourierTransformAtFrequency(k,w);

}

complex<double> CSCircleMaskGrating::surroundFourierTransformAtFrequency(vec2 k, double w) const
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


complex<double> CSCircleMaskGrating::centerFourierTransformAtFrequency(vec2 k, double w) const
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



double CSCircleMaskGrating::surroundMaskSize() const
{
    return m_surroundMaskSize;
}

double CSCircleMaskGrating::surroundOrientation(bool inDegrees) const
{
    if(inDegrees){
        return m_surroundOrientation * 180. / core::pi;
    }
    return m_surroundOrientation;

}

void CSCircleMaskGrating::setSurroundOrientation(double surroundOrientation)
{
    m_surroundOrientation = surroundOrientation;
}

double CSCircleMaskGrating::surroundPhase(bool inDegrees) const
{
    if(inDegrees){
        return m_surroundPhase * 180. / core::pi;
    }
    return m_surroundPhase;

}

double CSCircleMaskGrating::surroundContrast() const
{
    return m_surroundContrast;
}

double CSCircleMaskGrating::surroundSpatialFreq() const
{
    return m_surroundk;
}

void CSCircleMaskGrating::setSurroundSpatialFreq(double surroundSpatialFreq)
{
    m_surroundk = surroundSpatialFreq;
}

vec2 CSCircleMaskGrating::surroundkVec() const
{
    return m_surroundkVec;
}


double CSCircleMaskGrating::maskSize() const
{
    return m_maskSize;
}

double CSCircleMaskGrating::surroundTemporalFreq() const
{
    return m_surroundw;
}

