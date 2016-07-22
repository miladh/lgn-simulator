#ifndef CSCIRCLEMASKGRATING_H
#define CSCIRCLEMASKGRATING_H

#include "grating.h"

namespace lgnSimulator {
class CSCircleMaskGrating : public Grating
{
public:
    CSCircleMaskGrating(Integrator* const integrator,
                        double contrast, double phase,
                        double orientation, double maskSize,
                        double spatialFreq, double temporalFreq,
                        double surroundSpatialFreq, double surroundTemporalFreq,
                        double surroundContrast, double surroundPhase,
                        double surroundOrientation, double surroundMaskSize);
    ~CSCircleMaskGrating();

    // Stimulus interface
//    virtual void computeFourierTransform();

    void setSurroundSpatialFreq(double surroundSpatialFreq);
    void setSurroundOrientation(double surroundOrientation);

    double surroundSpatialFreq() const;
    double surroundContrast() const;
    double surroundPhase(bool inDegrees=false) const;
    double surroundOrientation(bool inDegrees=false) const;
    double surroundMaskSize() const;
    double maskSize() const;

    vec2 surroundkVec() const;

    double surroundTemporalFreq() const;


private:
    virtual double valueAtPoint(vec2 rVec, double t) const override;
    virtual complex<double> fourierTransformAtFrequency(vec2 k, double w) const override;
    complex<double> centerFourierTransformAtFrequency(vec2 k, double w) const;
    complex<double> surroundFourierTransformAtFrequency(vec2 k, double w) const;


    double m_surroundk = 0.0;
    double m_surroundw = 0.0;
    double m_surroundContrast = 1.0;
    double m_surroundPhase = 0;
    double m_surroundOrientation = 0;
    double m_surroundMaskSize = 0.0;
    double m_maskSize = 0.0;

    vec2 m_surroundkVec={0,0};

};
}

#endif // CSCIRCLEMASKGRATING_H
