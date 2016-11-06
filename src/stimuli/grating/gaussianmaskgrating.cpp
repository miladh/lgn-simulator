#include "gaussianmaskgrating.h"

/*!
 \class lgnSimulator::GaussianMaskGrating
 \inmodule lgnSimulator
 \ingroup lgnSimulator-gratingStimulus
 \brief Grating stimulus with gaussian mask.
 */

using namespace lgnSimulator;
GaussianMaskGrating::GaussianMaskGrating(Integrator* const integrator,
                                         double spatialFreq, double temporalFreq,
                                         double contrast, double phase,
                                         double orientation, double maskSize)
    : Grating(integrator, spatialFreq, temporalFreq, contrast, phase, orientation)
{
    m_mask= "gaussian";
    m_gaussianMask = new SpatialGaussian(maskSize); //WARNING! pointer is deleted.
}

GaussianMaskGrating::~GaussianMaskGrating()
{

}

double GaussianMaskGrating::valueAtPoint(vec2 r, double t) const
{

    double s = m_contrast * cos(dot(m_kVec, r) - m_w * t + m_phase)
            * m_gaussianMask->spatial(r);

    return s;

}

complex<double> GaussianMaskGrating::fourierTransformAtFrequency(vec2 k, double w) const
{

    complex<double> term1 = Special::delta(w, m_w) * m_phase_p_ft
            * m_gaussianMask->fourierTransform(m_kVec - k);

    complex<double> term2 = Special::delta(w, -m_w) * m_phase_m_ft
            * m_gaussianMask->fourierTransform(-m_kVec - k);

    return  core::pi * m_contrast * (term1 + term2) / m_integrator->temporalFreqResolution();

}
