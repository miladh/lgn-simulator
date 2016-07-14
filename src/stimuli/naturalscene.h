#ifndef NATURALSCENE_H
#define NATURALSCENE_H

#include "stimuli/stimulus.h"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/photo/photo.hpp"

namespace lgnSimulator {
class NaturalScene : public Stimulus
{
public:
    NaturalScene(const Integrator &integrator, string scenePath);
    ~NaturalScene();


    // Stimulus interface
public:
    void computeSpatiotemporal();
    void computeFourierTransform();


private:
    string m_sceneFilename;
    mat m_scene;
    cx_mat m_sceneFourierTransform;

    void readScene();

    virtual double temporalValueAtPoint(double t) = 0;
    virtual double fourierTransformAtTemporalFrequency(double w) = 0;
};
}
#endif // NATURALSCENE_H
