#ifndef CIRCLEMASKGRATING_H
#define CIRCLEMASKGRATING_H

#include "grating.h"

namespace lgnSimulator {
class CircleMaskGrating : public Grating
{
public:
    CircleMaskGrating(const Integrator &integrator,
                 double spatialFreq, double orientation, double temporalFreq, double contrast, double maskSize);
    ~CircleMaskGrating();

        // Grating interface
private:
    virtual double valueAtPoint(vec2 rVec, double t) const;
    virtual complex<double> fourierTransformAtFrequency(vec2 kVec, double w) const;

};

}
#endif // CIRCLEMASKGRATING_H
