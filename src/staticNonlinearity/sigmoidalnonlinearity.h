#ifndef SIGMOIDALNONLINEARITY_H
#define SIGMOIDALNONLINEARITY_H

#include "staticnonlinearity.h"
namespace lgnSimulator {
class SigmoidalNonlinearity : public StaticNonlinearity
{
public:
    SigmoidalNonlinearity(double maxValue, double halfMaxValue, double exponent);
    ~SigmoidalNonlinearity();

    // StaticNonlinearity interface
protected:
    virtual double advance(const double u) const;

private:
    double m_maxValue = 0.0;
    double m_halfmaxValue = 0.0;
    double m_exponent = 0.0;
};

}
lgnSimulator::SigmoidalNonlinearity createSigmoidalNonlinearity(const YAML::Node &cfg);

#endif // SIGMOIDALNONLINEARITY_H
