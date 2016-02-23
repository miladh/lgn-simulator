#include "special.h"

using namespace lgnSimulator;


Special::Special()
{

}

Special::~Special()
{

}

double Special::heaviside(double x)
{
    if (x < 0){
        return 0;
    }else{
        return 1;
    }
}

double Special::secondKindBesselFunction(double x)
{
    double j =  boost::math::cyl_bessel_j(1, x);
    return j;
}


double Special::delta(double x, double y) {
    if(x == y){
        return 1;
    }else{
        return 0;
    }
}

int Special::isOdd(int num)
{
    if(num % 2){
        return 1;
    }else{
        return 0;
    }
}



