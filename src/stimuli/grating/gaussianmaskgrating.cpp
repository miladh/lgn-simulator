#include "gaussianmaskgrating.h"

using namespace edog;


GaussianMaskGrating::GaussianMaskGrating(Integrator *integrator,
                                         vec2 kd, double wd, double contrast, double maskSize)

    : Grating(integrator, kd, wd, contrast, maskSize)
    , m_dog(new DOG(1,maskSize, 0, 1))
{
    m_mask = "gauss";
}

GaussianMaskGrating::~GaussianMaskGrating()
{

}


double GaussianMaskGrating::valueAtPoint(vec2 rVec, double t)
{

    double s = m_contrast * cos(dot(m_k, rVec) - m_w * t) * m_dog->spatial(rVec);
    return s;
}

double GaussianMaskGrating::fourierTransformAtFrequency(vec2 k, double w)
{
    double s = m_dog->fourierTransform(k)
            * Functions::delta(k[0], m_k[0])
            * Functions::delta(k[1], m_k[1])
            * Functions::delta(-w, m_w)
            /m_integrator->spatialFreqResolution()
            /m_integrator->spatialFreqResolution()
            /m_integrator->temporalFreqResolution();

    return 8*PI*PI*PI* m_contrast * s;
}


