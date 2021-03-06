#ifndef THRESHOLDNONLINEARITY_H
#define THRESHOLDNONLINEARITY_H

#include "staticnonlinearity.h"
#include "helper/special.h"


namespace lgnSimulator {
class ThresholdNonlinearity : public StaticNonlinearity
{
public:
    ThresholdNonlinearity(double threshold, double weight);
    ~ThresholdNonlinearity();

    // StaticNonlinearity interface
private:
    virtual double advance(const double u) const;
    double m_threshold = 0.0;
    double m_weight = 0.0;

};
}

lgnSimulator::ThresholdNonlinearity createThresholdNonlinearity(const YAML::Node &cfg);

#endif // THRESHOLDNONLINEARITY_H
