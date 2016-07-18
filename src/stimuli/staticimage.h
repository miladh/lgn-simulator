#ifndef STATICIMAGE_H
#define STATICIMAGE_H

#include "stimuli/naturalscene.h"


namespace lgnSimulator {
class StaticImage : public NaturalScene
{
public:
    StaticImage(Integrator* const integrator, string sceneFilename);
    ~StaticImage();

    // NaturalScene interface
private:
    virtual double temporalValueAtPoint(double t);
    virtual double fourierTransformAtTemporalFrequency(double w);
};

}
unique_ptr<lgnSimulator::StaticImage> createStaticImageStimulus(const lgnSimulator::Integrator &integrator, const YAML::Node &cfg);

#endif // STATICIMAGE_H
