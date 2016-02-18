#ifndef NATURALSCENEVIDEO_H
#define NATURALSCENEVIDEO_H

#include "stimuli/stimuli.h"

#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/photo/photo.hpp"
namespace lgnSimulator {
class NaturalSceneVideo : public Stimulus
{
public:
    NaturalSceneVideo(const Integrator &integrator, string sceneFilename);
    ~NaturalSceneVideo();

    // Stimulus interface
public:
    void computeSpatiotemporal();
    void computeFourierTransform();

private:
    string m_sceneFilename;
    cx_cube m_scene;

    void readScene();
};

}
unique_ptr<lgnSimulator::NaturalSceneVideo> createNaturalSceneVideoStimulus(const lgnSimulator::Integrator &integrator,
                                                  const YAML::Node &cfg);
#endif // NATURALSCENEVIDEO_H
