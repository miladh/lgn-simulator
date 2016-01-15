#ifndef STATICIMAGE_H
#define STATICIMAGE_H

#include "stimuli/naturalscene.h"


namespace edog {
class StaticImage : public NaturalScene
{
public:
    StaticImage(Integrator *integrator, string sceneFilename);
    ~StaticImage();

    // NaturalScene interface
private:
    virtual double temporalValueAtPoint(double t);
    virtual double fourierTransformAtTemporalFrequency(double w);
};

}
unique_ptr<edog::StaticImage> createStaticImageStimulus(edog::Integrator *integrator, const YAML::Node *cfg);

#endif // STATICIMAGE_H
