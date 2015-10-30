#ifndef STATICIMAGE_H
#define STATICIMAGE_H

#include "naturalscene.h"



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


StaticImage createStaticImageStimulus(Integrator *integrator, const Config *cfg);
#endif // STATICIMAGE_H
