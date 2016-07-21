#ifndef FULLFIELDGRATING_H
#define FULLFIELDGRATING_H

#include "grating.h"


namespace lgnSimulator {
class FullFieldGrating : public Grating
{
public:
    FullFieldGrating(Integrator* const integrator,
                     double spatialFreq, double temporalFreq,
                     double contrast, double phase,
                     double orientation);
    ~FullFieldGrating();

    // Grating interface
private:
    virtual double valueAtPoint(vec2 r, double t) const override;
    virtual complex<double> fourierTransformAtFrequency(vec2 k, double w) const override;

    double m_peak = 1.;
};

}
#endif // FULLFIELDGRATING_H
