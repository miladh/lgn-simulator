#ifndef NATURALSCENE_H
#define NATURALSCENE_H

#include "stimuli.h"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/photo/photo.hpp"


class NaturalScene : public Stimulus
{
public:
    NaturalScene(Integrator *integrator, string scenePath);
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
    string m_sceneFilename;
    cx_mat m_scene;
    cx_mat m_sceneFourierTransform;

    void readScene();

    virtual double temporalValueAtPoint(double t) = 0;
    virtual double fourierTransformAtTemporalFrequency(double w) = 0;
};

#endif // NATURALSCENE_H
