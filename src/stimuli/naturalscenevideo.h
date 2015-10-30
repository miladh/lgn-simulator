#ifndef NATURALSCENEVIDEO_H
#define NATURALSCENEVIDEO_H

#include "stimuli.h"

#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/photo/photo.hpp"

class NaturalSceneVideo : public Stimulus
{
public:
    NaturalSceneVideo(Integrator *integrator, string sceneFilename);
    ~NaturalSceneVideo();

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
    cx_cube m_scene;

    void readScene();
};

NaturalSceneVideo createNaturalSceneVideoStimulus(Integrator *integrator,
                                                  const Config *cfg);

#endif // NATURALSCENEVIDEO_H
