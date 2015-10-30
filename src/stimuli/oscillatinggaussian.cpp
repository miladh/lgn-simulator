#include "oscillatinggaussian.h"

OscillatingGaussian::OscillatingGaussian(Integrator *integrator,
                                         double exponent, double wd)
    : Stimulus(integrator)
    , m_exponent(exponent)
    , m_wd(wd)
{

}

OscillatingGaussian::~OscillatingGaussian()
{

}


double OscillatingGaussian::valueAtPoint(vec2 rVec, double t)
{
    return exp(-m_exponent*dot(rVec,rVec))*cos(m_wd * t);

}

double OscillatingGaussian::fourierTransformAtFrequency(vec2 k, double w)
{
   return PI/m_exponent*exp(-dot(k,k)/4./m_exponent)
           * Functions::delta(m_wd, w)*2*PI
           /m_integrator->temporalFreqResolution();
}

void OscillatingGaussian::computeSpatiotemporal()
{
    computeSpatiotemporalAnalytic();
}

void OscillatingGaussian::computeFourierTransform()
{
    computeFourierTransformAnalytic();
}


OscillatingGaussian createOscillatingGaussianStimulus(Integrator *integrator,
                                                      const Config *cfg)
{
    const Setting & root = cfg->getRoot();
    double exponent = root["stimuliSettings"]["exponent"];

    vec w = integrator->temporalFreqVec();
    double wd = w(1);

    return OscillatingGaussian(integrator, exponent, wd);
}

