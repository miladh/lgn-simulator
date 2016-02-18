#ifndef HEAVISIDENONLINEARITY_H
#define HEAVISIDENONLINEARITY_H

#include "staticnonlinearity.h"
#include "helper/specialfunctions.h"

namespace lgnSimulator {
class HeavisideNonlinearity : public StaticNonlinearity
{
public:
    HeavisideNonlinearity();
    ~HeavisideNonlinearity();

    // StaticNonlinearity interface
protected:
    virtual double advance(const double u) const;
};
}
#endif // HEAVISIDENONLINEARITY_H
