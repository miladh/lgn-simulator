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

SigmoidalNonlinearity createSigmoidalNonlinearity(const Config *cfg)
{
    const Setting & root = cfg->getRoot();
    double maxValue = root["staticNonlinearitySettings"]["maxValue"];
    double halfMaxValue = root["staticNonlinearitySettings"]["halfMaxValue"];
    double exponent = root["staticNonlinearitySettings"]["exponent"];

    return SigmoidalNonlinearity(maxValue, halfMaxValue, exponent);
}
