#ifndef SPECIALFUNCTIONS_H
#define SPECIALFUNCTIONS_H


#include <boost/math/special_functions/bessel.hpp>
#include <armadillo>

using namespace arma;

namespace lgnSimulator {
class SpecialFunctions
{
public:
    SpecialFunctions();
    ~SpecialFunctions();
    static double heaviside(double x);
    static double secondKindBesselFunction(double x);
    static double delta(double x, double y);
    static int isOdd(int num);

};
}
#endif // SPECIALFUNCTIONS_H
