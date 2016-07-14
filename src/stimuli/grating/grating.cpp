#include "grating.h"

#include "fullfieldgrating.h"
#include "circlemaskgrating.h"
#include "gaussianmaskgrating.h"
#include "cscirclemaskgrating.h"

using namespace lgnSimulator;


Grating::Grating(const Integrator &integrator,
                 double spatialFreq, double orientation, double temporalFreq,
                 double contrast, double maskSize)
    : Stimulus(integrator)
    , m_k(spatialFreq)
    , m_orientation(orientation*core::pi/180.)
    , m_w(temporalFreq)
    , m_contrast(contrast)
    , m_maskSize(maskSize)
{
    m_type = "grating";
    m_kVec = {Special::nearestValue(m_spatialFreqs, m_k*cos(m_orientation)),
              Special::nearestValue(m_spatialFreqs,  m_k*sin(m_orientation))};

    setSpatialFreq(sqrt(dot(m_kVec, m_kVec)));
    setOrientation(atan2(m_kVec(1), m_kVec(0)));

    if(maskSize > m_spatialVec.max()-m_spatialVec.min()){
        cerr << "Warning: mask size (" << maskSize
             << ") larger than grid length: " << m_spatialVec.max()-m_spatialVec.min()
             << endl;
    }
}

Grating::~Grating()
{
}


void Grating::computeSpatiotemporal()
{
    for(int k = 0; k < int(m_spatiotemporal.n_slices); k++){
        for(int i = 0; i < int(m_spatiotemporal.n_rows); i++){
            for(int j = 0; j < int(m_spatiotemporal.n_cols); j++){
                m_spatiotemporal(i,j,k) = valueAtPoint({m_spatialVec[i],
                                                        m_spatialVec[j]},
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


void Grating::setSpatialFreq(double spatialFreq)
{
    m_k = spatialFreq;
}

void Grating::setOrientation(double orientation)
{
    m_orientation = orientation;
}



double Grating::temporalFreq() const
{
    return m_w;
}

vec2 Grating::kVec() const
{
    return m_kVec;
}

string Grating::mask() const
{
    return m_mask;
}

double Grating::spatialFreq() const
{
    return m_k;
}

double Grating::orientation(bool inDegrees) const
{
    if(inDegrees){
        return m_orientation * 180. / core::pi;
    }
    return m_orientation;
}

double Grating::contrast() const
{
    return m_contrast;
}

double Grating::maskSize() const
{
    return m_maskSize;
}



unique_ptr<Grating> createGratingStimulus(const Integrator &integrator, const YAML::Node &cfg)
{
    string mask = cfg["mask"].as<string>();
    double maskSize = cfg["maskSize"].as<double>();
    double orientation = cfg["orientation"].as<double>();
    double contrast = cfg["C"].as<double>();
    int kId = cfg["kId"].as<int>();
    int wId = cfg["wId"].as<int>();

    vec k = integrator.spatialFreqVec();
    vec w = integrator.temporalFreqVec();


    if((kId  < -int(k.n_elem)/2) || (kId  > int(k.n_elem)/2-1)){
        cerr << "Too high or low index, kxId: " << kId << endl
             << "kId range: [" << -int(k.n_elem)/2 << "," << k.n_elem/2-1 <<"]" <<endl;
        return 0;

    }if(kId  < 0){
        kId+= k.n_elem;
    }

    if((wId  < -int(w.n_elem)/2) || (wId  > int(w.n_elem)/2-1)){
        cerr << "Too high or low index, wId: " << wId << endl
             << "wId range: [" << -int(w.n_elem)/2 << "," << w.n_elem/2-1 <<"]" <<endl;
        return 0;

    }if(wId  < 0){
        wId+= w.n_elem;
    }

    double temporalFreq = -w(wId); // -1 factor due to form of w vector
    double spatialFreq = k(kId);



    cout << "Stimulus: Grating with " << mask << " mask" << endl
         << "w=" << temporalFreq << endl;

    if(mask == "none"){
        return unique_ptr<FullFieldGrating>(
                    new FullFieldGrating(integrator, spatialFreq, orientation,
                                         temporalFreq, contrast));

    }else if(mask == "circle"){
        return unique_ptr<CircleMaskGrating>(
                    new CircleMaskGrating(integrator, spatialFreq, orientation,
                                          temporalFreq, contrast, maskSize));

    }else if(mask == "gaussian"){
        return unique_ptr<GaussianMaskGrating>(
                    new GaussianMaskGrating(integrator, spatialFreq, orientation,
                                          temporalFreq, contrast, maskSize));
    }else if(mask == "cscircle"){
        return unique_ptr<CSCircleMaskGrating>(
                    new CSCircleMaskGrating(integrator, spatialFreq, orientation,
                                          temporalFreq, contrast, maskSize));
    }else{
        cout << "mask: " << mask << endl;
        throw overflow_error("Unknown grating mask");
    }

}
