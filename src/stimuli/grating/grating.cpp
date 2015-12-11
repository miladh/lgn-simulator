#include "grating.h"

#include "fullfieldgrating.h"
#include "circlemaskgrating.h"
#include "gaussianmaskgrating.h"

Grating::Grating(Integrator *integrator,
                 vec2 kd, double wd, double contrast, string mask, double maskSize)
    : Stimulus(integrator)
    , m_k(kd)
    , m_w(wd)
    , m_contrast(contrast)
    , m_maskSize(maskSize * m_coordinateVec.max() * 2)
{

}

Grating::~Grating()
{
}


void Grating::computeSpatiotemporal()
{
    computeSpatiotemporalAnalytic();

}

void Grating::computeFourierTransform()
{
    computeFourierTransformAnalytic();
}


Grating* setGratingType(Integrator *integrator,
                                vec2 kd, double wd, double contrast,
                                string mask, double maskSize)
{
    if(mask == "none"){
        return new FullFieldGrating(integrator, kd, wd, contrast);

    }else if(mask == "gauss"){
        return new GaussianMaskGrating(integrator, kd, wd, contrast, maskSize);

    }else if(mask == "circle"){
        return new CircleMaskGrating(integrator, kd, wd, contrast, maskSize);

    }else{
        cout << "mask: " << mask << endl;
        throw overflow_error("Unknown grating mask");
    }
}


Grating* createGratingStimulus(Integrator *integrator, const Config *cfg)
{
    const Setting & root = cfg->getRoot();
    string mask = root["stimuliSettings"]["GratingSettings"]["mask"];
    double maskSize = root["stimuliSettings"]["GratingSettings"]["maskSize"];
    double contrast = root["stimuliSettings"]["GratingSettings"]["C"];

    vec k = integrator->spatialFreqVec();
    vec w = integrator->temporalFreqVec();
    double wd = w(1);
    double kx = k(4);
    double ky = k(0);

    return setGratingType(integrator, {kx, ky}, wd, contrast, mask, maskSize);
}
