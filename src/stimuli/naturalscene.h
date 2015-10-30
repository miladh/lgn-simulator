#ifndef NATURALSCENE_H
#define NATURALSCENE_H

#include "stimuli.h"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/photo/photo.hpp"


class NaturalScene : public Stimulus
{
public:
    NaturalScene(Integrator *integrator, mat scene);
    ~NaturalScene();


    // Stimulus interface
public:
    void computeSpatiotemporal();
    void computeFourierTransform();

private:
    double valueAtPoint(vec2 rVec, double t){
        (void)rVec;
        (void)t;
    }
    double fourierTransformAtFrequency(vec2 kVec, double w){
        (void)kVec;
        (void)w;
    }

private:
    cx_mat m_scene;
    cx_mat m_sceneFourierTransform;
};

NaturalScene createNaturalSceneStimulus(Integrator *integrator, const Config *cfg);

#endif // NATURALSCENE_H
