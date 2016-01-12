#ifndef THRESHOLDNONLINEARITY_H
#define THRESHOLDNONLINEARITY_H

#include "staticnonlinearity.h"
#include "../math/functions.h"

class ThresholdNonlinearity : public StaticNonlinearity
{
public:
    ThresholdNonlinearity(double threshold, double weight);
    ~ThresholdNonlinearity();

    // StaticNonlinearity interface
private:
    virtual double advance(const double u);

private:
    double m_threshold = 0.0;
    double m_weight = 0.0;

};

ThresholdNonlinearity createThresholdNonlinearity(const Config *cfg);

#endif // THRESHOLDNONLINEARITY_H
