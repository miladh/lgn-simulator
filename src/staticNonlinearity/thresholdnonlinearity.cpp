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


ThresholdNonlinearity createThresholdNonlinearity(const Config *cfg)
{
    const Setting & root = cfg->getRoot();
    double threshold = root["staticNonlinearitySettings"]["threshold"];
    double weight = root["staticNonlinearitySettings"]["weight"];

    return ThresholdNonlinearity(threshold, weight);
}
