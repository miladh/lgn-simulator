#ifndef SPECIAL_H
#define SPECIAL_H


#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/sinc.hpp>
#include <boost/tr1/cmath.hpp>
#include <helper/helperconstants.h>
#include <armadillo>
#include <gsl/gsl_sf_hyperg.h>

using namespace arma;

namespace lgnSimulator {
class Special
{
public:
    Special();
    ~Special();
    static double heaviside(double x);
    static int factorial(int n);
    static double secondKindBessel(double x);
    static double sinc(double x);
    static double confluentHypergeometric(double a, double b, double x);
    static double rect(double t, double T);
    static int isOdd(int num);
    static double delta(double x, double y);
    static double delta(vec2 x, vec2 y);
    static double nearestValue(const vec x, const double value);

};
}
#endif // SPECIAL_H
