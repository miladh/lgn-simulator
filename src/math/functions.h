#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#define PI 3.141592653589793



#include <boost/math/special_functions/bessel.hpp>
#include <armadillo>


using namespace arma;

namespace lgnSimulator {
class Functions
{
public:
    Functions();
    ~Functions();
    static double heaviside(double x);
    static double secondKindBesselFunction(double x);
    static double delta(double x, double y);
    static int isOdd(int num);

};
}
#endif // FUNCTIONS_H
