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
    static cx_cube fftShift(cx_cube c);


    static vec fftshift(vec x);
    static mat fftshift(mat x);



};

#endif // FUNCTIONS_H
