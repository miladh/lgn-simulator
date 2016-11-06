#include "grating.h"

#include "fullfieldgrating.h"
#include "circlemaskgrating.h"
#include "gaussianmaskgrating.h"
#include "cscirclemaskgrating.h"

/*!
 \class lgnSimulator::Grating
 \inmodule lgnSimulator-stimulus
 \brief Virtual class for grating stimulus.
 */


using namespace lgnSimulator;


Grating::Grating(Integrator* const integrator,
                 double spatialFreq, double temporalFreq,
                 double contrast, double phase,
                 double orientation)
    : Stimulus(integrator)
    , m_k(spatialFreq)
    , m_w(temporalFreq)
    , m_contrast(contrast)
    , m_phase(phase*core::pi/180.)
    , m_orientation(orientation*core::pi/180.)
    , m_phase_p_ft(exp(core::i * m_phase))
    , m_phase_m_ft(exp(-core::i * m_phase))
{
    m_type = "grating";
    m_kVec = {Special::nearestValue(m_spatialFreqs, m_k*cos(m_orientation)),
              Special::nearestValue(m_spatialFreqs,  m_k*sin(m_orientation))};

    setSpatialFreq(sqrt(dot(m_kVec, m_kVec)));
    setOrientation(atan2(m_kVec(1), m_kVec(0)));

//    cout << "K=" << m_k << endl;

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


double Grating::phase(bool inDegrees) const
{
    if(inDegrees){
        return m_phase * 180. / core::pi;
    }
    return m_phase;
}

double Grating::contrast() const
{
    return m_contrast;
}




unique_ptr<Grating> createGratingStimulus(Integrator* const integrator, const YAML::Node &cfg)
{
    string mask = cfg["mask"].as<string>();
    double orientation = cfg["orientation"].as<double>();
    double phase = cfg["phase"].as<double>();
    double contrast = cfg["C"].as<double>();
    int kId = cfg["kId"].as<int>();
    int wId = cfg["wId"].as<int>();

    vec k = integrator->spatialFreqVec();
    vec w = integrator->temporalFreqVec();


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
                    new FullFieldGrating(integrator, spatialFreq, temporalFreq,
                                         contrast, phase, orientation));

    }else if(mask == "circle"){
        double maskSize = cfg["maskSize"].as<double>();
        return unique_ptr<CircleMaskGrating>(
                    new CircleMaskGrating(integrator, spatialFreq, temporalFreq,
                                          contrast, phase, orientation, maskSize));

    }else if(mask == "gaussian"){
        double maskSize = cfg["maskSize"].as<double>();
        return unique_ptr<GaussianMaskGrating>(
                    new GaussianMaskGrating(integrator, spatialFreq, temporalFreq,
                                            contrast, phase, orientation, maskSize));
    }else if(mask == "cscircle"){
        double maskSize = cfg["maskSize"].as<double>();
        int surroundkId = cfg["surroundkId"].as<int>();
        int surroundwId = cfg["surroundwId"].as<int>();
        double surroundContrast = cfg["surroundC"].as<double>();
        double surroundPhase = cfg["surroundPhase"].as<double>();
        double surroundOrientation = cfg["surroundOrientation"].as<double>();
        double surroundMaskSize = cfg["surroundMaskSize"].as<double>();


        if((surroundkId  < -int(k.n_elem)/2) || (surroundkId  > int(k.n_elem)/2-1)){
            cerr << "Too high or low index, kxId: " << surroundkId << endl
                 << "kId range: [" << -int(k.n_elem)/2 << "," << k.n_elem/2-1 <<"]" <<endl;
            return 0;

        }if(surroundkId  < 0){
            surroundkId+= k.n_elem;
        }

        if((surroundwId  < -int(w.n_elem)/2) || (surroundwId  > int(w.n_elem)/2-1)){
            cerr << "Too high or low index, wId: " << surroundwId << endl
                 << "wId range: [" << -int(w.n_elem)/2 << "," << w.n_elem/2-1 <<"]" <<endl;
            return 0;

        }if(surroundwId  < 0){
            surroundwId+= w.n_elem;
        }

        double surroundTemporalFreq = -w(surroundwId); // -1 factor due to form of w vector
        double surroundSpatialFreq = k(surroundkId);




        return unique_ptr<CSCircleMaskGrating>(
                    new CSCircleMaskGrating(integrator, spatialFreq, temporalFreq,
                                            contrast, phase, orientation, maskSize,
                                            surroundSpatialFreq,  surroundTemporalFreq,
                                            surroundContrast,surroundPhase,
                                            surroundOrientation,  surroundMaskSize));
    }else{
        cout << "mask: " << mask << endl;
        throw overflow_error("Unknown grating mask");
    }

}
