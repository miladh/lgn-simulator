#ifndef HEAVISIDENONLINEARITY_H
#define HEAVISIDENONLINEARITY_H

#include "staticnonlinearity.h"
#include "math/functions.h"

class HeavisideNonlinearity : public StaticNonlinearity
{
public:
    HeavisideNonlinearity();
    ~HeavisideNonlinearity();

    // StaticNonlinearity interface
protected:
    virtual double advance(const double u);
};

#endif // HEAVISIDENONLINEARITY_H
