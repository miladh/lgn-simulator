#ifndef STATICIMAGE_H
#define STATICIMAGE_H

#include "stimuli/naturalscene.h"


namespace lgnSimulator {
class StaticImage : public NaturalScene
{
public:
    StaticImage(Integrator* const integrator, string sceneFilename, double delay, double period);
    ~StaticImage();

    // NaturalScene interface
private:
    virtual double temporalValueAtPoint(double t);
    virtual complex<double> fourierTransformAtTemporalFrequency(double w);
    double m_delay=0.0;
    double m_period=0.0;
};

}
unique_ptr<lgnSimulator::StaticImage> createStaticImageStimulus(lgnSimulator::Integrator* const integrator,
                                                                const YAML::Node &cfg);

#endif // STATICIMAGE_H
