#ifndef NATURALSCENE_H
#define NATURALSCENE_H

#include "stimuli.h"
#include "opencv2/highgui/highgui.hpp"



class NaturalScene : public Stimulus
{
public:
    NaturalScene(Integrator *integrator, mat scene);
    ~NaturalScene();

public:
    void computeSpatiotemporal();
    void computeFourierTransform();

    // Stimulus interface
private:
    virtual double valueAtPoint(vec2 rVec, double t);
    virtual double fourierTransformAtFrequency(vec2 k, double w);

    cx_mat m_scene;
    cx_mat m_sceneFourierTransform;
};

NaturalScene createNaturalSceneStimulus(Integrator *integrator, const Config *cfg);

#endif // NATURALSCENE_H
