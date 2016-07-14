#ifndef CSCIRCLEMASKGRATING_H
#define CSCIRCLEMASKGRATING_H

#include "grating.h"

namespace lgnSimulator {
class CSCircleMaskGrating : public Grating
{
public:
    CSCircleMaskGrating(const Integrator &integrator,
                        double spatialFreq, double orientation, double temporalFreq,
                        double contrast, double surroundSize);
    ~CSCircleMaskGrating();


    virtual void computeFourierTransform();
    // Grating interface
private:
    virtual double valueAtPoint(vec2 rVec, double t) const override;
    virtual complex<double> fourierTransformAtFrequency(vec2 k, double w) const override;
    double m_centerSize =0.0;
    vec2 m_centerK={0,0};
};
}

#endif // CSCIRCLEMASKGRATING_H
