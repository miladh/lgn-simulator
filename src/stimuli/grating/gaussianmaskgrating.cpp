#include "gaussianmaskgrating.h"

using namespace lgnSimulator;


GaussianMaskGrating::GaussianMaskGrating(const Integrator &integrator,
                                         vec2 kd, double wd, double contrast, double maskSize)

    : Grating(integrator, kd, wd, contrast, maskSize)
    , m_dog(new DOG(1,maskSize, 0, 1))
{
    m_mask = "gauss";
}

GaussianMaskGrating::~GaussianMaskGrating()
{

}


double GaussianMaskGrating::valueAtPoint(vec2 rVec, double t) const
{

    double s = m_contrast * cos(dot(m_k, rVec) - m_w * t) * m_dog->spatial(rVec);
    return s;
}

complex<double> GaussianMaskGrating::fourierTransformAtFrequency(vec2 k, double w) const
{
    complex<double> s = m_dog->fourierTransform(k)
            * Special::delta(k[0], m_k[0])
            * Special::delta(k[1], m_k[1])
            * Special::delta(-w, m_w)
            /m_integrator.spatialFreqResolution()
            /m_integrator.spatialFreqResolution()
            /m_integrator.temporalFreqResolution();

    return 8*core::pi*core::pi*core::pi* m_contrast * s;
}


