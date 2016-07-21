#include "fullfieldgrating.h"

using namespace lgnSimulator;


FullFieldGrating::FullFieldGrating(Integrator* const integrator,
                                   double spatialFreq, double temporalFreq,
                                   double contrast, double phase,
                                   double orientation)
    : Grating(integrator, spatialFreq, temporalFreq, contrast, phase, orientation)
{
    m_mask = "none";
    m_peak = 1.0
            / m_integrator->temporalFreqResolution()
            / m_integrator->spatialFreqResolution()
            / m_integrator->spatialFreqResolution();

}

FullFieldGrating::~FullFieldGrating()
{
}

double FullFieldGrating::valueAtPoint(vec2 r, double t) const
{
    double s = m_contrast * cos(dot(m_kVec, r) - m_w * t + m_phase);
    return s;
}

complex<double> FullFieldGrating::fourierTransformAtFrequency(vec2 k, double w) const
{
    complex<double> s = Special::delta(k, m_kVec) * Special::delta(w, m_w) * exp(core::i * m_phase)
                + Special::delta(k, -m_kVec) * Special::delta(w, -m_w)* exp(-core::i * m_phase);

    return 4.*core::pi*core::pi*core::pi * m_contrast * m_peak * s;
}


