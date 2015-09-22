#include "functions.h"

Functions::Functions()
{

}

Functions::~Functions()
{

}


vec Functions::fftshift(vec x)
{
    vec shifted;
    for(int i = 0; i < int(shifted.n_rows); i++){
        shifted(i)  = x(i) * pow(-1, i);
    }
    return shifted;
}

mat Functions::fftshift(mat x)
{
    vec shifted;
    for(int i = 0; i < int(shifted.n_rows); i++){
        for(int j = 0; j < int(shifted.n_cols); j++){

            shifted(i,j)  = x(i,j) * pow(-1, i+j);
        }
    }
    return shifted;
}

cx_cube Functions::fftShift(cx_cube c)
{

    cx_cube mShifted = c;
    for(int k = 0; k < int(c.n_slices); k++){
        for(int i = 0; i < int(c.n_rows); i++){
            for(int j = 0; j < int(c.n_cols); j++){

                mShifted(i,j,k)  = c(i,j,k) * pow(-1, i+j+k);
            }
        }
    }

    return mShifted;
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



