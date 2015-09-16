#include "functions.h"

Functions::Functions()
{

}

Functions::~Functions()
{

}

cx_mat Functions::fftShift(cx_mat m)
{

    cx_mat mShifted = m;
    for(int i = 0; i < int(m.n_cols); i++){
        for(int j = 0; j < int(m.n_rows); j++){

            mShifted(i,j)  = m(i,j) * pow(-1, i+j);
        }
    }

    return mShifted;

}

mat Functions::fftShift(mat m)
{

    mat mShifted = m;
    for(int i = 0; i < int(m.n_cols); i++){
        for(int j = 0; j < int(m.n_rows); j++){

            mShifted(i,j)  = m(i,j) * pow(-1, i+j);
        }
    }

    return mShifted;

}

cx_cube Functions::fftShift3d(cx_cube c)
{

    cx_cube mShifted = c;
    for(int k = 0; k < int(c.n_slices); k++){
        for(int i = 0; i < int(c.n_cols); i++){
            for(int j = 0; j < int(c.n_rows); j++){

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

