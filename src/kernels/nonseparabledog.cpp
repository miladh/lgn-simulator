#include "nonseparabledog.h"

using namespace lgnSimulator;


NonseparableDOG::NonseparableDOG(double weight,
                                 double a, double b, double c,
                                 double cenLatencyAlpha, double cenLatencyBeta,
                                 double surLatencyAlpha, double surLatencyBeta,
                                 double delay)
    : Kernel(weight)
    , m_relativeStrength(c)
{

    //!! create function points at deleted memory!!!
    m_spatialCenter = new SpatialGaussian(a);
    m_spatialSurround = new SpatialGaussian(b);
    m_temporalCenter = new DOE (cenLatencyAlpha, cenLatencyBeta, delay );
    m_temporalSurround = new DOE(surLatencyAlpha, surLatencyBeta, delay);

}


double lgnSimulator::NonseparableDOG::spatiotemporal(vec2 r, double t) const
{
     return   m_weight
             * (m_spatialCenter->spatial(r)
             * m_temporalCenter->temporal(t)
             - m_relativeStrength
             * m_spatialSurround->spatial(r)
             * m_temporalSurround->temporal(t));
}

complex<double> lgnSimulator::NonseparableDOG::fourierTransform(vec2 k, double w) const
{
 return m_weight *
         (m_spatialCenter->fourierTransform(k)
          * m_temporalCenter->fourierTransform(w)
          - m_relativeStrength
          * m_spatialSurround->fourierTransform(k)
          * m_temporalSurround->fourierTransform(w));
}


NonseparableDOG createNonseparableDOGKernel(const YAML::Node &cfg)
{

    double weight = cfg["weight"].as<double>();

    double cenLatencyAlpha = cfg["cenLatencyAlpha"].as<double>();
    double cenLatencyBeta  = cfg["cenLatencyBeta"].as<double>();
    double surLatencyAlpha = cfg["surLatencyAlpha"].as<double>();
    double surLatencyBeta = cfg["surLatencyBeta"].as<double>();
    double delay    = cfg["delay"].as<double>();

    double a = cfg["a"].as<double>();
    double b = cfg["b"].as<double>();
    double c = cfg["c"].as<double>();


    return NonseparableDOG(weight, a, b, c,
                           cenLatencyAlpha, cenLatencyBeta,
                           surLatencyAlpha, surLatencyBeta,
                           delay);
}
