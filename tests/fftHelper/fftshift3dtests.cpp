#include <unittest++/UnitTest++.h>
#include <armadillo>

#include "math/ffthelper.h"


using namespace std;
using namespace arma;

SUITE(fftHelper){
    TEST(fftshift_3D){
        cx_cube sig, centered, shifted;


        sig.set_size( 4, 2, 3 );
        centered.set_size( 4, 2, 3 );

        sig.slice( 0 ) << -3 << 4 << endr
                       << 3 << -4 << endr
                       << -2 << 5 << endr
                       << -3 << 1 << endr;
        sig.slice( 1 ) << -5 << -3 << endr
                       << 5 << 2 << endr
                       << -5 << -1 << endr
                       << 2 << -2 << endr;
        sig.slice( 2 ) << 2 << 3 << endr
                       << 3 << 0 << endr
                       << -1 << 1 << endr
                       << 5 << -2 << endr;

        centered.slice( 0 ) << 1 << -1 << endr
                            << -2 << 5 << endr
                            << 3 << 2 << endr
                            << 0 << 3 << endr;
        centered.slice( 1 ) << 5 << -2 << endr
                            << 1 << -3 << endr
                            << 4 << -3 << endr
                            << -4 << 3 << endr;
        centered.slice( 2 ) << -1 << -5 << endr
                            << -2 << 2 << endr
                            << -3 << -5 << endr
                            << 2 << 5 << endr;

        shifted = FFTHelper::fftShift(sig);

        for(int k = 0; k < int(sig.n_slices); k++){
            for(int i = 0; i < int(sig.n_rows); i++){
                for(int j = 0; j < int(sig.n_cols); j++){
                    CHECK_EQUAL(real(shifted(i,j,k)), real(centered(i,j,k)));
                }
            }
        }
        /*------------------------------------------------------*/







    }


}
