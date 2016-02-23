#include "fullfieldgrating.h"

using namespace lgnSimulator;


FullFieldGrating::FullFieldGrating(const Integrator &integrator,
                                   vec2 kd, double wd, double contrast)
    : Grating(integrator, kd, wd, contrast, 0)
{
   m_mask = "none";
}

FullFieldGrating::~FullFieldGrating()
{

}


double FullFieldGrating::valueAtPoint(vec2 rVec, double t) const
{
    double s = m_contrast * cos(dot(m_k, rVec) - m_w * t);
    return s;
}

complex<double> FullFieldGrating::fourierTransformAtFrequency(vec2 k, double w) const
{
    double s = Special::delta(k[0], m_k[0])
            * Special::delta(k[1], m_k[1])
            * Special::delta(-w, m_w)
            /m_integrator.spatialFreqResolution()
            /m_integrator.spatialFreqResolution()
            /m_integrator.temporalFreqResolution();

    return 8*core::pi*core::pi*core::pi* m_contrast * s;
}
