#include "functions.h"

using namespace edog;


Functions::Functions()
{

}

Functions::~Functions()
{

}

double Functions::heaviside(double x)
{
    if (x < 0){
        return 0;
    }else{
        return 1;
    }
}

double Functions::secondKindBesselFunction(double x)
{
    double j =  boost::math::cyl_bessel_j(1, x);
    return j;
}


double Functions::delta(double x, double y) {
    if(x == y){
        return 1;
    }else{
        return 0;
    }
}

int Functions::isOdd(int num)
{
    if(num % 2){
        return 1;
    }else{
        return 0;
    }
}



