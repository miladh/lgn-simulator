#ifndef CIRCLEMASKGRATING_H
#define CIRCLEMASKGRATING_H

#include "grating.h"

namespace lgnSimulator {
class CircleMaskGrating : public Grating
{
public:
    CircleMaskGrating(Integrator* const integrator,
                      double spatialFreq, double temporalFreq,
                      double contrast, double phase,
                      double orientation, double maskSize);
    ~CircleMaskGrating();

    double maskSize() const;
    virtual complex<double> fourierTransformAtFrequency(vec2 k, double w) const override;

private:
    virtual double valueAtPoint(vec2 rVec, double t) const override;

    double m_maskSize = 0.0;
};

}
#endif // CIRCLEMASKGRATING_H
