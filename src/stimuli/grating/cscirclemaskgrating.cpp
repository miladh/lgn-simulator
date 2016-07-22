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
    , m_phase_ps_ft(exp(core::i * m_surroundPhase))
    , m_phase_ms_ft(exp(-core::i * m_surroundPhase))
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

    if(m_surroundMaskSize < m_maskSize){
        throw invalid_argument("m_surroundMaskSize < m_maskSize: "
                               +to_string(m_surroundMaskSize)+ "< "
                               + to_string(m_maskSize));
    }

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

    return (centerFourierTransformAtFrequency(k,w) +surroundFourierTransformAtFrequency(k,w))
            * core::pi * core::pi * 0.25 / m_integrator->temporalFreqResolution();

}

complex<double> CSCircleMaskGrating::surroundFourierTransformAtFrequency(vec2 k, double w) const
{
    complex<double> wm = Special::delta(w, m_surroundw) * m_phase_ps_ft;
    complex<double> wp = Special::delta(w, -m_surroundw)* m_phase_ms_ft;

    double km = sqrt(dot(k - m_surroundkVec, k - m_surroundkVec));
    double kp = sqrt(dot(k + m_surroundkVec, k + m_surroundkVec));

    double arg_km_ds = km * m_surroundMaskSize * 0.5;
    double arg_km_dc = km * m_maskSize * 0.5;
    double arg_kp_ds = kp * m_surroundMaskSize * 0.5;
    double arg_kp_dc = kp * m_maskSize * 0.5;


    double S_km_ds = m_surroundMaskSize * m_surroundMaskSize;
    double S_kp_ds = S_km_ds;
    double S_km_dc = m_maskSize * m_maskSize;
    double S_kp_dc = S_km_dc;

    if(arg_km_ds!= 0){
        S_km_ds *= 2 * Special::secondKindBessel(arg_km_ds)/arg_km_ds;
    }
    if(arg_kp_ds!= 0){
        S_kp_ds *= 2 * Special::secondKindBessel(arg_kp_ds)/arg_kp_ds;
    }
    if(arg_km_dc!= 0){
        S_km_dc *= 2 * Special::secondKindBessel(arg_km_dc)/arg_km_dc;
    }
    if(arg_kp_dc!= 0){
        S_kp_dc *= 2 * Special::secondKindBessel(arg_kp_dc)/arg_kp_dc;
    }

    return m_surroundContrast * (wm * (S_km_ds - S_km_dc) + wp* (S_kp_ds - S_kp_dc));

}


complex<double> CSCircleMaskGrating::centerFourierTransformAtFrequency(vec2 k, double w) const
{

    double arg1 = sqrt(dot(k - m_kVec, k - m_kVec))* m_maskSize * 0.5;
    double arg2 = sqrt(dot(k + m_kVec, k + m_kVec))* m_maskSize * 0.5;

    complex<double> term1 = Special::delta(w, m_w) * m_phase_p_ft;// CHECK SIGN!!!
    complex<double> term2 = Special::delta(w, -m_w) * m_phase_m_ft;

    if(arg1!= 0){
        term1 *= 2*Special::secondKindBessel(arg1)/arg1;
    }
    if(arg2!= 0){
        term2 *= 2*Special::secondKindBessel(arg2)/arg2;
    }

    return m_contrast * m_maskSize * m_maskSize * (term1 + term2);
}


//void CSCircleMaskGrating::computeFourierTransform()
//{
//    computeSpatiotemporal();
//    m_fourierTransform = m_integrator->forwardFFT(m_spatiotemporal);

//}



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

