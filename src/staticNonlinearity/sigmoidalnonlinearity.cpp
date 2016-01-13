#include "sigmoidalnonlinearity.h"

SigmoidalNonlinearity::SigmoidalNonlinearity(double maxValue, double halfMaxValue, double exponent)
    : m_maxValue(maxValue)
    , m_halfmaxValue(halfMaxValue)
    , m_exponent(exponent)
{
}

SigmoidalNonlinearity::~SigmoidalNonlinearity()
{

}



double SigmoidalNonlinearity::advance(const double u)
{
    return m_maxValue / (1 + exp(m_exponent * (u - m_halfmaxValue)));
}

SigmoidalNonlinearity createSigmoidalNonlinearity(const YAML::Node *cfg)
{

    double maxValue = (*cfg)["staticNonlinearitySettings"]["maxValue"].as<double>();
    double halfMaxValue = (*cfg)["staticNonlinearitySettings"]["halfMaxValue"].as<double>();
    double exponent = (*cfg)["staticNonlinearitySettings"]["exponent"].as<double>();

    return SigmoidalNonlinearity(maxValue, halfMaxValue, exponent);
}
