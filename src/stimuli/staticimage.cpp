#include "staticimage.h"

using namespace lgnSimulator;


StaticImage::StaticImage(Integrator* const integrator, string sceneFilename):
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
    return Special::delta(0, w)/m_integrator->temporalFreqResolution();
}


unique_ptr<StaticImage> createStaticImageStimulus(Integrator* const integrator, const YAML::Node &cfg)
{
    //Read file
    string sceneFilename = cfg["sceneFilename"].as<string>();

    return unique_ptr<StaticImage>(new StaticImage (integrator, sceneFilename));
}
