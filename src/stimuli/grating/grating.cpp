#include "grating.h"

#include "fullfieldgrating.h"
#include "circlemaskgrating.h"
#include "gaussianmaskgrating.h"

using namespace edog;


Grating::Grating(Integrator *integrator,
                 vec2 kd, double wd, double contrast, double maskSize)
    : Stimulus(integrator)
    , m_k(kd)
    , m_w(wd)
    , m_contrast(contrast)
    , m_maskSize(maskSize * m_coordinateVec.max() * 2)
{
    m_type = "grating";
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

double Grating::w() const
{
    return m_w;
}

vec2 Grating::k() const
{
    return m_k;
}



unique_ptr<Grating> createGratingStimulus(Integrator *integrator, const YAML::Node *cfg)
{
    string mask = (*cfg)["mask"].as<string>();
    double maskSize = (*cfg)["maskSize"].as<double>();
    double contrast = (*cfg)["C"].as<double>();
    int kxId = (*cfg)["kxId"].as<int>();
    int kyId = (*cfg)["kyId"].as<int>();
    int wId = (*cfg)["wId"].as<int>();

    vec k = integrator->spatialFreqVec();
    vec w = integrator->temporalFreqVec();


    if((kxId  < -int(k.n_elem)/2) || (kxId  > int(k.n_elem)/2-1)){
        cerr << "Too high index, kxId: " << kxId << endl
             << "kxId range: [" << -int(k.n_elem)/2 << "," << k.n_elem/2-1 <<"]" <<endl;
        return 0;

    }if(kxId  < 0){
        kxId+= k.n_elem;
    }

    if((kyId  < -int(k.n_elem)/2) || (kyId  > int(k.n_elem)/2-1)){
        cerr << "Too high index, kyId: " << kyId << endl
             << "kyId range: [" << -int(k.n_elem)/2 << "," << k.n_elem/2-1 <<"]" <<endl;
        return 0;

    }if(kyId  < 0){
        kyId+= k.n_elem;
    }

    if((wId  < -int(w.n_elem)/2) || (wId  > int(w.n_elem)/2-1)){
        cerr << "Too high index, wId: " << wId << endl
             << "wId range: [" << -int(w.n_elem)/2 << "," << w.n_elem/2-1 <<"]" <<endl;
        return 0;

    }if(wId  < 0){
        wId+= w.n_elem;
    }

    double wd = w(wId);
    double kx = k(kxId);
    double ky = k(kyId);

    cout << "Stimulus: Grating with " << mask << " mask" << endl
         << "kx=" << kx << ", ky=" << ky << ", w=" << wd << endl;



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
