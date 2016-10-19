#include "staticimage.h"

using namespace lgnSimulator;


StaticImage::StaticImage(Integrator* const integrator, string sceneFilename, double delay, double period)
    : NaturalScene(integrator, sceneFilename)
    , m_delay(delay)
    , m_period(period)
{
}

StaticImage::~StaticImage()
{

}

double StaticImage::temporalValueAtPoint(double t)
{
    return Special::heaviside(t-m_delay) /*Special::rect(t-m_delay, m_period)*/;
}

complex<double> StaticImage::fourierTransformAtTemporalFrequency(double w)
{

    double rePart=core::pi*Special::delta(0., w)/m_integrator->temporalFreqResolution();
    complex<double> imPart=1.0/(-core::i);
    if(w==0){
        imPart/=m_integrator->temporalFreqResolution();
    }else{
        imPart/=w;
    }
    return exp(core::i*w*m_delay)*(imPart + rePart);

//    return m_period * Special::sinc(w*m_period/2.)*exp(core::i*w*(m_delay+m_period/2));
}


unique_ptr<StaticImage> createStaticImageStimulus(Integrator* const integrator, const YAML::Node &cfg)
{
    //Read file
    string sceneFilename = cfg["sceneFilename"].as<string>();
    double delay = cfg["delay"].as<double>();
    double period = cfg["period"].as<double>();

    return unique_ptr<StaticImage>(new StaticImage (integrator, sceneFilename, delay, period));
}
