#include "thresholdnonlinearity.h"

ThresholdNonlinearity::ThresholdNonlinearity(double threshold, double weight)
    : m_threshold(threshold)
    , m_weight(weight)
{

}

ThresholdNonlinearity::~ThresholdNonlinearity()
{

}


double ThresholdNonlinearity::advance(const double u)
{
    double x = u - m_threshold;
    if(Functions::heaviside(x)){
        return x * m_weight;
    }else{
        return 0;
    }
}


ThresholdNonlinearity createThresholdNonlinearity(const YAML::Node *cfg)
{
    double threshold = (*cfg)["staticNonlinearitySettings"]["threshold"].as<double>();
    double weight = (*cfg)["staticNonlinearitySettings"]["weight"].as<double>();

    return ThresholdNonlinearity(threshold, weight);
}
