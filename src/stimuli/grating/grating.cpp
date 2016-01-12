#include "grating.h"

#include "fullfieldgrating.h"
#include "circlemaskgrating.h"
#include "gaussianmaskgrating.h"

Grating::Grating(Integrator *integrator,
                 vec2 kd, double wd, double contrast, double maskSize)
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
    for(int k = 0; k < int(m_spatioTemporal.n_slices); k++){
        for(int i = 0; i < int(m_spatioTemporal.n_rows); i++){
            for(int j = 0; j < int(m_spatioTemporal.n_cols); j++){
                m_spatioTemporal(i,j,k) = valueAtPoint({m_coordinateVec[i],
                                                        m_coordinateVec[j]},
                                                       m_timeVec[k]);
            }
        }
    }

}

void Grating::computeFourierTransform()
{
    for(int k = 0; k < int(m_fourierTransform.n_slices); k++){
        for(int i = 0; i < int(m_fourierTransform.n_rows); i++){
            for(int j = 0; j < int(m_fourierTransform.n_cols); j++){
                m_fourierTransform(i,j,k) =
                        fourierTransformAtFrequency({m_spatialFreqs[i],
                                                     m_spatialFreqs[j]},
                                                    m_temporalFreqs[k]);

            }
        }
    }
}



unique_ptr<Grating> createGratingStimulus(Integrator *integrator, const Config *cfg)
{
    const Setting & root = cfg->getRoot();
    string mask = root["stimuliSettings"]["GratingSettings"]["mask"];
    double maskSize = root["stimuliSettings"]["GratingSettings"]["maskSize"];
    double contrast = root["stimuliSettings"]["GratingSettings"]["C"];

    vec k = integrator->spatialFreqVec();
    vec w = integrator->temporalFreqVec();
    double wd = w(1);
    double kx = k(10);
    double ky = k(0);

    if(mask == "none"){
        return unique_ptr<FullFieldGrating>(new FullFieldGrating(integrator, {kx, ky}, wd, contrast));

    }else if(mask == "gauss"){
        return unique_ptr<GaussianMaskGrating>(new GaussianMaskGrating(integrator, {kx, ky}, wd, contrast, maskSize));

    }else if(mask == "circle"){
        return unique_ptr<CircleMaskGrating>(new CircleMaskGrating(integrator, {kx, ky}, wd, contrast, maskSize));

    }else{
        cout << "mask: " << mask << endl;
        throw overflow_error("Unknown grating mask");
    }
}
