#include "patchgrating.h"

PatchGrating::PatchGrating(Integrator *integrator, vec2 kd,
                           double wd, double contrast,
                           double spotDiameterRatio)
    : Stimulus(integrator)
    , m_k(kd)
    , m_w(wd)
    , m_contrast(contrast)
    , m_spotDiameter(spotDiameterRatio * m_coordinateVec.max() * 2)
{
}

PatchGrating::~PatchGrating()
{

}

double PatchGrating::valueAtPoint(vec2 rVec, double t)
{
    double r = sqrt(dot(rVec, rVec));
    double s = m_contrast * (1 - Functions::heaviside(r - m_spotDiameter * 0.5))
            * cos(dot(m_k, rVec) - m_w * t);

    return s;

}



double PatchGrating::fourierTransformAtFrequency(vec2 kVec, double w)
{
    if(!Functions::delta(w, -m_w)){
        return 0;
    }

    double s = m_contrast * PI * PI * m_spotDiameter * m_spotDiameter * 0.5;
    double arg = sqrt(dot(kVec - m_k, kVec - m_k)) * m_spotDiameter * 0.5;

    if(arg != 0){
        s *= 2. * Functions::secondKindBesselFunction(arg)/arg;
    }

    return s/m_integrator->temporalFreqResolution();
}



void PatchGrating::computeSpatiotemporal()
{
    computeSpatiotemporalAnalytic();

}

void PatchGrating::computeFourierTransform()
{
    computeFourierTransformAnalytic();
}



PatchGrating createPatchGratingStimulus(Integrator *integrator, const Config *cfg)
{
    const Setting & root = cfg->getRoot();
    double contrast = root["stimuliSettings"]["C"];
    double spotDiameterRatio = root["stimuliSettings"]["d"];

    vec k = integrator->spatialFreqVec();
    vec w = integrator->temporalFreqVec();
    double wd = w(1);
    double kx = k(4);
    double ky = k(0);


    return PatchGrating(integrator, {kx, ky}, wd, contrast, spotDiameterRatio);
}
