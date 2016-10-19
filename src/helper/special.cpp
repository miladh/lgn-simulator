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

int Special::factorial(int n)
{
    if(n > 12){
        cerr << "too large (n_max=12) n in factorial(n): " << n << endl;
    }
    return boost::math::factorial<double>(n);
}


double Special::confluentHypergeometric(double a, double b, double x)
{
    double F = gsl_sf_hyperg_1F1(a, b, x);
    return F;
}

double Special::rect(double t, double T)
{
    if(t <= T){
        return 1;
    }else{
        return 0;
    }

}



double Special::secondKindBessel(double x)
{
    double j = boost::math::cyl_bessel_j(1, x);
    return j;
}

double Special::sinc(double x)
{
    return  boost::math::sinc_pi(x);
}


double Special::delta(double x, double y) {
    if(fabs(x - y) < core::epsilon){
        return 1;
    }else{
        return 0;
    }
}

double Special::delta(vec2 x, vec2 y) {
       return delta(x(0), y(0)) * delta(x(1), y(1));
}

int Special::isOdd(int num)
{
    if(num % 2){
        return 1;
    }else{
        return 0;
    }
}


double Special::nearestValue(const vec x, const double value) {

    double nearestValue=INFINITY;
    double difference = INFINITY;

    for(double xi: x){
        double dx = fabs(xi - value);
        if(fabs(xi - value) < difference){
            difference = dx;
            nearestValue = xi;
        }
    }
       return nearestValue;
}


