#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <boost/math/special_functions/bessel.hpp>

class Functions
{
public:
    Functions();
    ~Functions();
    static double heaviside(double x);
    static double secondKindBesselFunction(double x);
    static double delta(double x, double y);
};

#endif // FUNCTIONS_H
