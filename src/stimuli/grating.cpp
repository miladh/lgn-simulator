#include "grating.h"

Grating::Grating(Integrator integrator, vec2 kd, double wd, double contrast)
    : Stimulus(integrator)
    , m_k(kd)
    , m_w(wd)
    , m_contrast(contrast)
{
}

Grating::~Grating()
{

}

double Grating::valueAtPoint(vec2 rVec, double t)
{

    double s = m_contrast * cos(dot(m_k, rVec) - m_w * t);
    return s;
}

double Grating::fourierTransformAtFrequency(vec2 k, double w)
{
    double s = Functions::delta(k[0], m_k[0])
            * Functions::delta(k[1], m_k[1])
            * Functions::delta(w, -m_w)
            /m_integrator.spatialFreqResolution()
            /m_integrator.spatialFreqResolution()
            /m_integrator.temporalFreqResolution();

    return 8*PI*PI*PI* m_contrast * s;
}


Grating *createGratingStimulus(Integrator integrator, const Config *cfg)
{
    const Setting & root = cfg->getRoot();
    double contrast = root["stimuliSettings"]["C"];

    vec k = integrator.spatialFreqVec();
    vec w = integrator.temporalFreqVec();
    double wd = w(w.n_elem/2+2);
    double kx = k(k.n_elem/2+6);
    double ky = k(k.n_elem/2);

    return new Grating(integrator, {kx, ky}, wd, contrast);
}
