#ifndef SPECIAL_H
#define SPECIAL_H


#include <boost/math/special_functions/bessel.hpp>
#include <armadillo>

using namespace arma;

namespace lgnSimulator {
class Special
{
public:
    Special();
    ~Special();
    static double heaviside(double x);
    static double secondKindBesselFunction(double x);
    static double delta(double x, double y);
    static int isOdd(int num);

};
}
#endif // SPECIAL_H
