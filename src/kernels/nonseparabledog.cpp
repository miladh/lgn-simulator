#include "nonseparabledog.h"

using namespace lgnSimulator;

NonseparableDOG::NonseparableDOG(double weight,
                                 double A, double a, double B, double b,
                                 double cenLatencyAlpha, double cenLatencyBeta,
                                 double surLatencyAlpha, double surLatencyBeta,
                                 double delay):
    Kernel(weight)
{
    m_spatialCentre = new SpatialGaussian(A,a);
    m_spatialSurround = new SpatialGaussian(B,b);
    m_temporalCenter = new DOE (cenLatencyAlpha, cenLatencyBeta, delay );
    m_temporalSurround = new DOE(surLatencyAlpha, surLatencyBeta, delay);

}


double lgnSimulator::NonseparableDOG::spatiotemporal(vec2 r, double t) const
{
     return m_spatialCentre->spatial(r) * m_temporalCenter->temporal(t)
             - m_spatialSurround->spatial(r) * m_temporalSurround->temporal(t);
}

complex<double> lgnSimulator::NonseparableDOG::fourierTransform(vec2 k, double w) const
{
 return m_spatialCentre->fourierTransform(k) * m_temporalCenter->fourierTransform(w)
         - m_spatialSurround->fourierTransform(k) * m_temporalSurround->fourierTransform(w);
}


NonseparableDOG createNonseparableDOGKernel(const YAML::Node &cfg)
{

    double weight = cfg["weight"].as<double>();

    double cenLatencyAlpha = cfg["cenLatencyAlpha"].as<double>();
    double cenLatencyBeta  = cfg["cenLatencyBeta"].as<double>();
    double surLatencyAlpha = cfg["surLatencyAlpha"].as<double>();
    double surLatencyBeta = cfg["surLatencyBeta"].as<double>();
    double delay    = cfg["delay"].as<double>();

    double A = cfg["A"].as<double>();
    double a = cfg["a"].as<double>();
    double B = cfg["B"].as<double>();
    double b = cfg["b"].as<double>();


    return NonseparableDOG(weight, A, a, B, b,
                           cenLatencyAlpha, cenLatencyBeta,
                           surLatencyAlpha, surLatencyBeta,
                           delay);
}
