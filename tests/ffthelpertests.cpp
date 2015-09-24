#include <unittest++/UnitTest++.h>
#include <armadillo>


#include "math/ffthelper.h"


using namespace std;
using namespace arma;


SUITE(fftHelper){
    TEST(fftshift_1D){

        cx_vec A = {0.,  1.,  2.,  3.,  4., -5., -4., -3., -2., -1.};
        cx_vec centered = {-5., -4., -3., -2., -1.,  0.,  1.,  2.,  3.,  4.};

        cx_vec B = FFTHelper::fftShift(A);
        for(int i = 0; i < int(A.n_elem); i++){
            CHECK_EQUAL(B(i), centered(i));
        }
    }

    TEST(fftshift_2D){

        cx_mat A;
        A << 0.0  << 1.0  << 2.0  << endr
          << 3.0  << 4.0  << -4.0 << endr
          << -3.0 << -2.0 << -1.0 << endr;

        cx_mat centered;
        centered << -1.0  << -3.0  << -2.0  << endr
                 << 2.0  << 0.0  << 1.0 << endr
                 << -4.0 << 3.0 << 4.0 << endr;

        cx_mat B = FFTHelper::fftShift(A);

        for(int i = 0; i < int(A.n_rows); i++){
            for(int j = 0; j < int(A.n_cols); j++){
                CHECK_EQUAL(B(i,j), centered(i,j));
            }
        }



    }

    TEST(fftshift_3D){

        cx_cube A(3,3,2);
        A.slice(0) << 0.0  << 1.0  << 2.0  << endr
                   << 3.0  << 4.0  << -4.0 << endr
                   << -3.0 << -2.0 << -1.0 << endr;

        A.slice(1) << 10.0  << 11.0  << 12.0  << endr
                   << 13.0  << 14.0  << -14.0 << endr
                   << -13.0 << -12.0 << -11.0 << endr;


        cx_cube centered= 0*A;
        centered.slice(0) << -11.0  << -13.0  << -12.0  << endr
                          << 12.0  << 10.0  << 11.0 << endr
                          << -14.0 << 13.0 << 14.0 << endr;

        centered.slice(1) << -1.0  << -3.0  << -2.0  << endr
                          << 2.0  << 0.0  << 1.0 << endr
                          << -4.0 << 3.0 << 4.0 << endr;



        cx_cube B = FFTHelper::fftShift(A);

        for(int k = 0; k < int(A.n_slices); k++){
            for(int i = 0; i < int(A.n_rows); i++){
                for(int j = 0; j < int(A.n_cols); j++){
                    CHECK_EQUAL(B(i,j,k), centered(i,j,k));
                }
            }
        }

    }


    TEST(fftshift_3D_1){

        cx_cube A(3,2,2);
        A.slice(0) << 0.0  << 1.0  << endr
                   << 3.0  << 4.0  << endr
                   << -3.0 << -2.0  << endr;

        A.slice(1) << 10.0  << 11.0   << endr
                   << 13.0  << 14.0   << endr
                   << -13.0 << -12.0  << endr;


        cx_cube centered= 0*A;
        centered.slice(0) << -12.0  << -13.0   << endr
                          << 11.0  << 10.0  << endr
                          << 14.0 << 13.0  << endr;

        centered.slice(1) << -2.0  << -3.0   << endr
                          << 1.0  << 0.0  <<  endr
                          << 4.0 << 3.0 << endr;


        cx_cube B = FFTHelper::fftShift(A);

        for(int k = 0; k < int(A.n_slices); k++){
            for(int i = 0; i < int(A.n_rows); i++){
                for(int j = 0; j < int(A.n_cols); j++){
                    CHECK_EQUAL(B(i,j,k), centered(i,j,k));
                }
            }
        }


    }
}
