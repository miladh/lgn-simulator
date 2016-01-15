#ifndef NATURALSCENEVIDEO_H
#define NATURALSCENEVIDEO_H

#include "stimuli/stimuli.h"

#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/photo/photo.hpp"
namespace edog {
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
    string m_sceneFilename;
    cx_cube m_scene;

    void readScene();
};

}
unique_ptr<edog::NaturalSceneVideo> createNaturalSceneVideoStimulus(edog::Integrator *integrator,
                                                  const YAML::Node *cfg);
#endif // NATURALSCENEVIDEO_H
