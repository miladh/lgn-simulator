#ifndef CIRCLEMASKGRATING_H
#define CIRCLEMASKGRATING_H

#include "grating.h"

namespace lgnSimulator {
class CircleMaskGrating : public Grating
{
public:
    CircleMaskGrating(Integrator* const integrator,
                      double spatialFreq, double orientation, double temporalFreq,
                      double contrast, double maskSize, double phase=0.0);
    ~CircleMaskGrating();

//    virtual void computeFourierTransform();
        // Grating interface
private:
    virtual double valueAtPoint(vec2 rVec, double t) const override;
    virtual complex<double> fourierTransformAtFrequency(vec2 k, double w) const override;
};

}
#endif // CIRCLEMASKGRATING_H
