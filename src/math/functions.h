#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#define PI 3.14159265359



#include <boost/math/special_functions/bessel.hpp>
#include <armadillo>


using namespace arma;

class Functions
{
public:
    Functions();
    ~Functions();
    static double heaviside(double x);
    static double secondKindBesselFunction(double x);
    static double delta(double x, double y);
    static cx_mat fftShift(cx_mat m);
    static mat fftShift(mat m);
    static cx_cube fftShift3d(cx_cube c);


};

#endif // FUNCTIONS_H
