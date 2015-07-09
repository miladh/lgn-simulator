#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <boost/math/special_functions/bessel.hpp>

class Functions
{
public:
    Functions();
    ~Functions();
    double heaviside(double x);
    double secondKindBesselFunction(double x);
    double delta(double x, double y);
};

#endif // FUNCTIONS_H
