#include "thresholdnonlinearity.h"

using namespace lgnSimulator;


ThresholdNonlinearity::ThresholdNonlinearity(double threshold, double weight)
    : m_threshold(threshold)
    , m_weight(weight)
{

}

ThresholdNonlinearity::~ThresholdNonlinearity()
{

}


double ThresholdNonlinearity::advance(const double u) const
{
    double x = u - m_threshold;
    if(SpecialFunctions::heaviside(x)){
        return x * m_weight;
    }else{
        return 0;
    }
}


ThresholdNonlinearity createThresholdNonlinearity(const YAML::Node &cfg)
{
    double threshold = cfg["threshold"].as<double>();
    double weight = cfg["weight"].as<double>();

    return ThresholdNonlinearity(threshold, weight);
}
