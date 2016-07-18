#include "gaussianmaskgrating.h"

using namespace lgnSimulator;
GaussianMaskGrating::GaussianMaskGrating(Integrator* const integrator,
                                         double spatialFreq, double orientation,
                                         double temporalFreq,double contrast,
                                         double maskSize)
    : Grating(integrator, spatialFreq, orientation, temporalFreq, contrast, maskSize)
{
    m_mask= "gaussian";
    m_gaussianMask = new SpatialGaussian(maskSize);
}

GaussianMaskGrating::~GaussianMaskGrating()
{

}

double GaussianMaskGrating::valueAtPoint(vec2 r, double t) const
{

    double s = m_contrast * cos(dot(m_kVec, r) - m_w * t)
            * m_gaussianMask->spatial(r);

    return s;

}

complex<double> GaussianMaskGrating::fourierTransformAtFrequency(vec2 k, double w) const
{

    complex<double> term1 = Special::delta(w, m_w)
            * m_gaussianMask->fourierTransform(m_kVec - k);

    complex<double> term2 = Special::delta(w, -m_w)
            * m_gaussianMask->fourierTransform(-m_kVec - k);

    return  core::pi * m_contrast * (term1 + term2) / m_integrator->temporalFreqResolution();

}
