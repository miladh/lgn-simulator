#include "staticimage.h"

StaticImage::StaticImage(Integrator *integrator, string sceneFilename):
    NaturalScene(integrator, sceneFilename)
{
}

StaticImage::~StaticImage()
{

}

double StaticImage::temporalValueAtPoint(double t)
{
    (void) t;
    return 1.;
}

double StaticImage::fourierTransformAtTemporalFrequency(double w)
{
    return Functions::delta(0, w)/m_integrator->temporalFreqResolution();
}


unique_ptr<StaticImage> createStaticImageStimulus(Integrator *integrator, const Config *cfg)
{
    //Read file
    const Setting & root = cfg->getRoot();
    string sceneFilename = root["stimuliSettings"]["sceneFilename"];

    return unique_ptr<StaticImage>(new StaticImage (integrator, sceneFilename));
}
