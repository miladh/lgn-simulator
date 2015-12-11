#include "gaussianmaskgrating.h"

GaussianMaskGrating::GaussianMaskGrating(Integrator *integrator,
                                         vec2 kd, double wd, double contrast, double maskSize)
    : Grating(integrator, kd, wd, contrast, maskSize)
{
}

GaussianMaskGrating::~GaussianMaskGrating()
{

}


double GaussianMaskGrating::valueAtPoint(vec2 rVec, double t)
{

    double s = m_contrast * cos(dot(m_k, rVec) - m_w * t);
    return s;
}

double GaussianMaskGrating::fourierTransformAtFrequency(vec2 k, double w)
{
    double s = Functions::delta(k[0], m_k[0])
            * Functions::delta(k[1], m_k[1])
            * Functions::delta(-w, m_w)
            /m_integrator->spatialFreqResolution()
            /m_integrator->spatialFreqResolution()
            /m_integrator->temporalFreqResolution();

    return 8*PI*PI*PI* m_contrast * s;
}




GaussianMaskGrating createGaussianMaskGratingStimulus(Integrator *integrator, const Config *cfg)
{
    const Setting & root = cfg->getRoot();
    double contrast = root["stimuliSettings"]["GratingSettings"]["C"];
    double maskSize = root["stimuliSettings"]["GratingSettings"]["maskSize"];

    vec k = integrator->spatialFreqVec();
    vec w = integrator->temporalFreqVec();
    double wd = w(1);
    double kx = k(4);
    double ky = k(0);

    return GaussianMaskGrating(integrator, {kx, ky}, wd, contrast, maskSize);
}
