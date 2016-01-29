#ifndef NATURALSCENE_H
#define NATURALSCENE_H

#include "stimuli/stimuli.h"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/photo/photo.hpp"

namespace edog {
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
    string m_sceneFilename;
    cx_mat m_scene;
    cx_mat m_sceneFourierTransform;

    void readScene();

    virtual double temporalValueAtPoint(double t) = 0;
    virtual double fourierTransformAtTemporalFrequency(double w) = 0;
};
}
#endif // NATURALSCENE_H