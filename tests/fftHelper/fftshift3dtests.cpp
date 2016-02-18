#include <unittest++/UnitTest++.h>
#include <armadillo>

#include "helper/ffthelper.h"


using namespace std;
using namespace arma;
using namespace lgnSimulator;

SUITE(fftHelper){
    TEST(fftshift_3D){
        cx_cube sig, centered, shifted;


        sig.set_size( 5, 3, 3 );
        centered.set_size( 5, 3, 3 );
        sig.slice( 0 ) << -1 << -2 << -4 << endr
                        << 0 << 1 << 4 << endr
                        << -4 << 2 << 1 << endr
                        << -3 << -5 << -3 << endr
                        << 4 << 3 << 0 << endr;
        sig.slice( 1 ) << 3 << -5 << 4 << endr
                        << 1 << 3 << 3 << endr
                        << -2 << -4 << 1 << endr
                        << 2 << 0 << 0 << endr
                        << -2 << 1 << 1 << endr;
        sig.slice( 2 ) << -2 << -1 << -2 << endr
                        << -2 << -1 << -4 << endr
                        << 2 << -3 << -1 << endr
                        << -5 << -4 << -5 << endr
                        << -4 << 4 << -3 << endr;

        centered.slice( 0 ) << -5 << -5 << -4 << endr
                              << -3 << -4 << 4 << endr
                              << -2 << -2 << -1 << endr
                              << -4 << -2 << -1 << endr
                              << -1 << 2 << -3 << endr;
        centered.slice( 1 ) << -3 << -3 << -5 << endr
                              << 0 << 4 << 3 << endr
                              << -4 << -1 << -2 << endr
                              << 4 << 0 << 1 << endr
                              << 1 << -4 << 2 << endr;
        centered.slice( 2 ) << 0 << 2 << 0 << endr
                              << 1 << -2 << 1 << endr
                              << 4 << 3 << -5 << endr
                              << 3 << 1 << 3 << endr
                              << 1 << -2 << -4 << endr;

        shifted = FFTHelper::fftShift(sig);
        for(int k = 0; k < int(sig.n_slices); k++){
            for(int i = 0; i < int(sig.n_rows); i++){
                for(int j = 0; j < int(sig.n_cols); j++){
                     CHECK_EQUAL(real(shifted(i,j,k)), real(centered(i,j,k)));
                }
             }
        }
        /*------------------------------------------------------*/

        sig.set_size( 3, 2, 5 );
        centered.set_size( 3, 2, 5 );
        sig.slice( 0 ) << -2 << 3 << endr
                        << 1 << -2 << endr
                        << 5 << -1 << endr;
        sig.slice( 1 ) << 4 << 2 << endr
                        << -1 << -4 << endr
                        << -3 << -1 << endr;
        sig.slice( 2 ) << 4 << -4 << endr
                        << 3 << -2 << endr
                        << -1 << -5 << endr;
        sig.slice( 3 ) << -4 << -5 << endr
                        << -1 << -5 << endr
                        << 3 << -3 << endr;
        sig.slice( 4 ) << -5 << -3 << endr
                        << 4 << 0 << endr
                        << -1 << -5 << endr;

        centered.slice( 0 ) << -3 << 3 << endr
                              << -5 << -4 << endr
                              << -5 << -1 << endr;
        centered.slice( 1 ) << -5 << -1 << endr
                              << -3 << -5 << endr
                              << 0 << 4 << endr;
        centered.slice( 2 ) << -1 << 5 << endr
                              << 3 << -2 << endr
                              << -2 << 1 << endr;
        centered.slice( 3 ) << -1 << -3 << endr
                              << 2 << 4 << endr
                              << -4 << -1 << endr;
        centered.slice( 4 ) << -5 << -1 << endr
                              << -4 << 4 << endr
                              << -2 << 3 << endr;

        shifted = FFTHelper::fftShift(sig);
        for(int k = 0; k < int(sig.n_slices); k++){
            for(int i = 0; i < int(sig.n_rows); i++){
                for(int j = 0; j < int(sig.n_cols); j++){
                     CHECK_EQUAL(real(shifted(i,j,k)), real(centered(i,j,k)));
                }
             }
        }
        /*------------------------------------------------------*/

        sig.set_size( 4, 3, 5 );
        centered.set_size( 4, 3, 5 );
        sig.slice( 0 ) << 4 << -5 << -5 << endr
                        << -4 << 5 << 5 << endr
                        << 1 << 3 << -4 << endr
                        << 3 << -3 << -4 << endr;
        sig.slice( 1 ) << 3 << 4 << 3 << endr
                        << 2 << 2 << -1 << endr
                        << 1 << -3 << 4 << endr
                        << -1 << -5 << 5 << endr;
        sig.slice( 2 ) << -1 << 5 << 4 << endr
                        << 0 << -3 << 3 << endr
                        << -2 << -5 << 2 << endr
                        << 2 << 0 << 5 << endr;
        sig.slice( 3 ) << -2 << -5 << -4 << endr
                        << 3 << -1 << -3 << endr
                        << -5 << 1 << 2 << endr
                        << -2 << -1 << -4 << endr;
        sig.slice( 4 ) << 0 << 2 << -1 << endr
                        << 5 << 4 << -3 << endr
                        << -5 << -3 << 1 << endr
                        << -5 << 2 << 5 << endr;

        centered.slice( 0 ) << 2 << -5 << 1 << endr
                              << -4 << -2 << -1 << endr
                              << -4 << -2 << -5 << endr
                              << -3 << 3 << -1 << endr;
        centered.slice( 1 ) << 1 << -5 << -3 << endr
                              << 5 << -5 << 2 << endr
                              << -1 << 0 << 2 << endr
                              << -3 << 5 << 4 << endr;
        centered.slice( 2 ) << -4 << 1 << 3 << endr
                              << -4 << 3 << -3 << endr
                              << -5 << 4 << -5 << endr
                              << 5 << -4 << 5 << endr;
        centered.slice( 3 ) << 4 << 1 << -3 << endr
                              << 5 << -1 << -5 << endr
                              << 3 << 3 << 4 << endr
                              << -1 << 2 << 2 << endr;
        centered.slice( 4 ) << 2 << -2 << -5 << endr
                              << 5 << 2 << 0 << endr
                              << 4 << -1 << 5 << endr
                              << 3 << 0 << -3 << endr;

        shifted = FFTHelper::fftShift(sig);
        for(int k = 0; k < int(sig.n_slices); k++){
            for(int i = 0; i < int(sig.n_rows); i++){
                for(int j = 0; j < int(sig.n_cols); j++){
                     CHECK_EQUAL(real(shifted(i,j,k)), real(centered(i,j,k)));
                }
             }
        }
        /*------------------------------------------------------*/

        sig.set_size( 4, 2, 2 );
        centered.set_size( 4, 2, 2 );
        sig.slice( 0 ) << -2 << 0 << endr
                        << -4 << 2 << endr
                        << -4 << 1 << endr
                        << -4 << 5 << endr;
        sig.slice( 1 ) << -1 << -4 << endr
                        << -2 << -3 << endr
                        << -1 << 3 << endr
                        << 2 << 5 << endr;

        centered.slice( 0 ) << 3 << -1 << endr
                              << 5 << 2 << endr
                              << -4 << -1 << endr
                              << -3 << -2 << endr;
        centered.slice( 1 ) << 1 << -4 << endr
                              << 5 << -4 << endr
                              << 0 << -2 << endr
                              << 2 << -4 << endr;

        shifted = FFTHelper::fftShift(sig);
        for(int k = 0; k < int(sig.n_slices); k++){
            for(int i = 0; i < int(sig.n_rows); i++){
                for(int j = 0; j < int(sig.n_cols); j++){
                     CHECK_EQUAL(real(shifted(i,j,k)), real(centered(i,j,k)));
                }
             }
        }
        /*------------------------------------------------------*/

        sig.set_size( 2, 3, 5 );
        centered.set_size( 2, 3, 5 );
        sig.slice( 0 ) << 3 << -2 << -1 << endr
                        << 4 << -2 << 4 << endr;
        sig.slice( 1 ) << 4 << -2 << -5 << endr
                        << 5 << 4 << 1 << endr;
        sig.slice( 2 ) << 2 << -5 << -5 << endr
                        << 3 << 2 << 3 << endr;
        sig.slice( 3 ) << -5 << 0 << 5 << endr
                        << 0 << -3 << 0 << endr;
        sig.slice( 4 ) << 2 << -4 << 4 << endr
                        << 4 << 4 << 1 << endr;

        centered.slice( 0 ) << 0 << 0 << -3 << endr
                              << 5 << -5 << 0 << endr;
        centered.slice( 1 ) << 1 << 4 << 4 << endr
                              << 4 << 2 << -4 << endr;
        centered.slice( 2 ) << 4 << 4 << -2 << endr
                              << -1 << 3 << -2 << endr;
        centered.slice( 3 ) << 1 << 5 << 4 << endr
                              << -5 << 4 << -2 << endr;
        centered.slice( 4 ) << 3 << 3 << 2 << endr
                              << -5 << 2 << -5 << endr;

        shifted = FFTHelper::fftShift(sig);
        for(int k = 0; k < int(sig.n_slices); k++){
            for(int i = 0; i < int(sig.n_rows); i++){
                for(int j = 0; j < int(sig.n_cols); j++){
                     CHECK_EQUAL(real(shifted(i,j,k)), real(centered(i,j,k)));
                }
             }
        }
        /*------------------------------------------------------*/

        sig.set_size( 3, 2, 2 );
        centered.set_size( 3, 2, 2 );
        sig.slice( 0 ) << 1 << 2 << endr
                        << -1 << -1 << endr
                        << 1 << 0 << endr;
        sig.slice( 1 ) << 0 << 2 << endr
                        << -1 << 3 << endr
                        << 1 << 5 << endr;

        centered.slice( 0 ) << 5 << 1 << endr
                              << 2 << 0 << endr
                              << 3 << -1 << endr;
        centered.slice( 1 ) << 0 << 1 << endr
                              << 2 << 1 << endr
                              << -1 << -1 << endr;

        shifted = FFTHelper::fftShift(sig);
        for(int k = 0; k < int(sig.n_slices); k++){
            for(int i = 0; i < int(sig.n_rows); i++){
                for(int j = 0; j < int(sig.n_cols); j++){
                     CHECK_EQUAL(real(shifted(i,j,k)), real(centered(i,j,k)));
                }
             }
        }
        /*------------------------------------------------------*/

        sig.set_size( 2, 5, 5 );
        centered.set_size( 2, 5, 5 );
        sig.slice( 0 ) << 5 << 2 << 1 << 1 << 2 << endr
                        << -2 << -5 << -1 << -2 << -4 << endr;
        sig.slice( 1 ) << -2 << -2 << 0 << 0 << 1 << endr
                        << 1 << 1 << 3 << 5 << -1 << endr;
        sig.slice( 2 ) << -2 << -5 << -2 << -1 << -2 << endr
                        << -2 << -2 << 3 << 2 << -4 << endr;
        sig.slice( 3 ) << -4 << 2 << -2 << 1 << 4 << endr
                        << -2 << 4 << -2 << 4 << 5 << endr;
        sig.slice( 4 ) << 1 << 4 << 4 << -5 << -1 << endr
                        << -3 << -1 << -2 << -4 << 5 << endr;

        centered.slice( 0 ) << 4 << 5 << -2 << 4 << -2 << endr
                              << 1 << 4 << -4 << 2 << -2 << endr;
        centered.slice( 1 ) << -4 << 5 << -3 << -1 << -2 << endr
                              << -5 << -1 << 1 << 4 << 4 << endr;
        centered.slice( 2 ) << -2 << -4 << -2 << -5 << -1 << endr
                              << 1 << 2 << 5 << 2 << 1 << endr;
        centered.slice( 3 ) << 5 << -1 << 1 << 1 << 3 << endr
                              << 0 << 1 << -2 << -2 << 0 << endr;
        centered.slice( 4 ) << 2 << -4 << -2 << -2 << 3 << endr
                              << -1 << -2 << -2 << -5 << -2 << endr;

        shifted = FFTHelper::fftShift(sig);
        for(int k = 0; k < int(sig.n_slices); k++){
            for(int i = 0; i < int(sig.n_rows); i++){
                for(int j = 0; j < int(sig.n_cols); j++){
                     CHECK_EQUAL(real(shifted(i,j,k)), real(centered(i,j,k)));
                }
             }
        }
        /*------------------------------------------------------*/

        sig.set_size( 4, 3, 5 );
        centered.set_size( 4, 3, 5 );
        sig.slice( 0 ) << 4 << 4 << -5 << endr
                        << 0 << -5 << -5 << endr
                        << 4 << 4 << 1 << endr
                        << -2 << 5 << 2 << endr;
        sig.slice( 1 ) << 0 << -5 << -4 << endr
                        << 4 << 5 << -5 << endr
                        << -1 << -5 << 2 << endr
                        << -5 << 5 << -2 << endr;
        sig.slice( 2 ) << 1 << 5 << 3 << endr
                        << -3 << -1 << -2 << endr
                        << -3 << -2 << -3 << endr
                        << 4 << 0 << 1 << endr;
        sig.slice( 3 ) << -4 << 2 << 5 << endr
                        << 1 << 0 << 1 << endr
                        << 3 << 1 << -3 << endr
                        << 5 << 5 << 3 << endr;
        sig.slice( 4 ) << 1 << -3 << -3 << endr
                        << 5 << -5 << 2 << endr
                        << 1 << 0 << 2 << endr
                        << 0 << 5 << 3 << endr;

        centered.slice( 0 ) << -3 << 3 << 1 << endr
                              << 3 << 5 << 5 << endr
                              << 5 << -4 << 2 << endr
                              << 1 << 1 << 0 << endr;
        centered.slice( 1 ) << 2 << 1 << 0 << endr
                              << 3 << 0 << 5 << endr
                              << -3 << 1 << -3 << endr
                              << 2 << 5 << -5 << endr;
        centered.slice( 2 ) << 1 << 4 << 4 << endr
                              << 2 << -2 << 5 << endr
                              << -5 << 4 << 4 << endr
                              << -5 << 0 << -5 << endr;
        centered.slice( 3 ) << 2 << -1 << -5 << endr
                              << -2 << -5 << 5 << endr
                              << -4 << 0 << -5 << endr
                              << -5 << 4 << 5 << endr;
        centered.slice( 4 ) << -3 << -3 << -2 << endr
                              << 1 << 4 << 0 << endr
                              << 3 << 1 << 5 << endr
                              << -2 << -3 << -1 << endr;

        shifted = FFTHelper::fftShift(sig);
        for(int k = 0; k < int(sig.n_slices); k++){
            for(int i = 0; i < int(sig.n_rows); i++){
                for(int j = 0; j < int(sig.n_cols); j++){
                     CHECK_EQUAL(real(shifted(i,j,k)), real(centered(i,j,k)));
                }
             }
        }
        /*------------------------------------------------------*/

        sig.set_size( 3, 3, 2 );
        centered.set_size( 3, 3, 2 );
        sig.slice( 0 ) << 4 << -1 << 2 << endr
                        << -4 << 1 << 3 << endr
                        << -2 << -5 << 5 << endr;
        sig.slice( 1 ) << -5 << 0 << 0 << endr
                        << 3 << 4 << 1 << endr
                        << 3 << 4 << -3 << endr;

        centered.slice( 0 ) << -3 << 3 << 4 << endr
                              << 0 << -5 << 0 << endr
                              << 1 << 3 << 4 << endr;
        centered.slice( 1 ) << 5 << -2 << -5 << endr
                              << 2 << 4 << -1 << endr
                              << 3 << -4 << 1 << endr;

        shifted = FFTHelper::fftShift(sig);
        for(int k = 0; k < int(sig.n_slices); k++){
            for(int i = 0; i < int(sig.n_rows); i++){
                for(int j = 0; j < int(sig.n_cols); j++){
                     CHECK_EQUAL(real(shifted(i,j,k)), real(centered(i,j,k)));
                }
             }
        }
        /*------------------------------------------------------*/

        sig.set_size( 4, 3, 2 );
        centered.set_size( 4, 3, 2 );
        sig.slice( 0 ) << 2 << 2 << 5 << endr
                        << 1 << -3 << 1 << endr
                        << 0 << 3 << 5 << endr
                        << -1 << 1 << -5 << endr;
        sig.slice( 1 ) << 2 << 3 << 3 << endr
                        << 2 << 0 << 1 << endr
                        << 5 << 2 << -1 << endr
                        << -4 << 0 << 1 << endr;

        centered.slice( 0 ) << -1 << 5 << 2 << endr
                              << 1 << -4 << 0 << endr
                              << 3 << 2 << 3 << endr
                              << 1 << 2 << 0 << endr;
        centered.slice( 1 ) << 5 << 0 << 3 << endr
                              << -5 << -1 << 1 << endr
                              << 5 << 2 << 2 << endr
                              << 1 << 1 << -3 << endr;

        shifted = FFTHelper::fftShift(sig);
        for(int k = 0; k < int(sig.n_slices); k++){
            for(int i = 0; i < int(sig.n_rows); i++){
                for(int j = 0; j < int(sig.n_cols); j++){
                     CHECK_EQUAL(real(shifted(i,j,k)), real(centered(i,j,k)));
                }
             }
        }
        /*------------------------------------------------------*/

        sig.set_size( 5, 5, 4 );
        centered.set_size( 5, 5, 4 );
        sig.slice( 0 ) << -2 << -1 << 2 << -1 << 5 << endr
                        << -5 << -4 << 2 << 3 << -2 << endr
                        << 4 << -2 << -5 << 2 << -5 << endr
                        << -4 << 4 << 5 << -4 << -3 << endr
                        << -5 << -3 << 5 << 5 << -3 << endr;
        sig.slice( 1 ) << -2 << 0 << -5 << 3 << -1 << endr
                        << -2 << -5 << -4 << 3 << -4 << endr
                        << 4 << -1 << -1 << 2 << -4 << endr
                        << -1 << 5 << -5 << 2 << -4 << endr
                        << -3 << 5 << 3 << 4 << -1 << endr;
        sig.slice( 2 ) << -3 << 4 << 5 << 4 << -3 << endr
                        << -1 << 5 << 2 << -2 << 3 << endr
                        << 5 << 0 << -1 << -1 << 5 << endr
                        << -4 << -1 << 5 << 4 << 5 << endr
                        << -4 << 1 << -4 << 5 << -1 << endr;
        sig.slice( 3 ) << 1 << 1 << 4 << 2 << -4 << endr
                        << -5 << -5 << -5 << 4 << -4 << endr
                        << -2 << 5 << 1 << 2 << -5 << endr
                        << -2 << -1 << -1 << -2 << -1 << endr
                        << -4 << 4 << -4 << 2 << 2 << endr;

        centered.slice( 0 ) << 4 << 5 << -4 << -1 << 5 << endr
                              << 5 << -1 << -4 << 1 << -4 << endr
                              << 4 << -3 << -3 << 4 << 5 << endr
                              << -2 << 3 << -1 << 5 << 2 << endr
                              << -1 << 5 << 5 << 0 << -1 << endr;
        centered.slice( 1 ) << -2 << -1 << -2 << -1 << -1 << endr
                              << 2 << 2 << -4 << 4 << -4 << endr
                              << 2 << -4 << 1 << 1 << 4 << endr
                              << 4 << -4 << -5 << -5 << -5 << endr
                              << 2 << -5 << -2 << 5 << 1 << endr;
        centered.slice( 2 ) << -4 << -3 << -4 << 4 << 5 << endr
                              << 5 << -3 << -5 << -3 << 5 << endr
                              << -1 << 5 << -2 << -1 << 2 << endr
                              << 3 << -2 << -5 << -4 << 2 << endr
                              << 2 << -5 << 4 << -2 << -5 << endr;
        centered.slice( 3 ) << 2 << -4 << -1 << 5 << -5 << endr
                              << 4 << -1 << -3 << 5 << 3 << endr
                              << 3 << -1 << -2 << 0 << -5 << endr
                              << 3 << -4 << -2 << -5 << -4 << endr
                              << 2 << -4 << 4 << -1 << -1 << endr;

        shifted = FFTHelper::fftShift(sig);
        for(int k = 0; k < int(sig.n_slices); k++){
            for(int i = 0; i < int(sig.n_rows); i++){
                for(int j = 0; j < int(sig.n_cols); j++){
                     CHECK_EQUAL(real(shifted(i,j,k)), real(centered(i,j,k)));
                }
             }
        }
        /*------------------------------------------------------*/

        sig.set_size( 4, 3, 5 );
        centered.set_size( 4, 3, 5 );
        sig.slice( 0 ) << -4 << -3 << 1 << endr
                        << 4 << 2 << -5 << endr
                        << -2 << -2 << 3 << endr
                        << -5 << 1 << 1 << endr;
        sig.slice( 1 ) << 3 << -5 << -1 << endr
                        << 0 << -3 << -4 << endr
                        << -4 << 0 << -4 << endr
                        << -2 << 0 << -5 << endr;
        sig.slice( 2 ) << 1 << -1 << 0 << endr
                        << -4 << -2 << 0 << endr
                        << 0 << -3 << 2 << endr
                        << -1 << -4 << 2 << endr;
        sig.slice( 3 ) << 0 << 1 << -1 << endr
                        << -5 << -3 << 4 << endr
                        << 1 << -3 << -1 << endr
                        << 1 << 3 << -3 << endr;
        sig.slice( 4 ) << -5 << 5 << -5 << endr
                        << 2 << 4 << 2 << endr
                        << 0 << 1 << -4 << endr
                        << -5 << -2 << -2 << endr;

        centered.slice( 0 ) << -1 << 1 << -3 << endr
                              << -3 << 1 << 3 << endr
                              << -1 << 0 << 1 << endr
                              << 4 << -5 << -3 << endr;
        centered.slice( 1 ) << -4 << 0 << 1 << endr
                              << -2 << -5 << -2 << endr
                              << -5 << -5 << 5 << endr
                              << 2 << 2 << 4 << endr;
        centered.slice( 2 ) << 3 << -2 << -2 << endr
                              << 1 << -5 << 1 << endr
                              << 1 << -4 << -3 << endr
                              << -5 << 4 << 2 << endr;
        centered.slice( 3 ) << -4 << -4 << 0 << endr
                              << -5 << -2 << 0 << endr
                              << -1 << 3 << -5 << endr
                              << -4 << 0 << -3 << endr;
        centered.slice( 4 ) << 2 << 0 << -3 << endr
                              << 2 << -1 << -4 << endr
                              << 0 << 1 << -1 << endr
                              << 0 << -4 << -2 << endr;

        shifted = FFTHelper::fftShift(sig);
        for(int k = 0; k < int(sig.n_slices); k++){
            for(int i = 0; i < int(sig.n_rows); i++){
                for(int j = 0; j < int(sig.n_cols); j++){
                     CHECK_EQUAL(real(shifted(i,j,k)), real(centered(i,j,k)));
                }
             }
        }
        /*------------------------------------------------------*/

        sig.set_size( 3, 4, 3 );
        centered.set_size( 3, 4, 3 );
        sig.slice( 0 ) << -5 << -1 << 3 << -3 << endr
                        << -1 << 3 << 2 << 1 << endr
                        << 0 << -2 << -2 << 3 << endr;
        sig.slice( 1 ) << -4 << -5 << -4 << -2 << endr
                        << 5 << 2 << -4 << -1 << endr
                        << 2 << -4 << 2 << 5 << endr;
        sig.slice( 2 ) << 0 << 4 << 5 << -5 << endr
                        << 1 << 1 << -1 << 5 << endr
                        << -2 << 4 << -3 << 3 << endr;

        centered.slice( 0 ) << -3 << 3 << -2 << 4 << endr
                              << 5 << -5 << 0 << 4 << endr
                              << -1 << 5 << 1 << 1 << endr;
        centered.slice( 1 ) << -2 << 3 << 0 << -2 << endr
                              << 3 << -3 << -5 << -1 << endr
                              << 2 << 1 << -1 << 3 << endr;
        centered.slice( 2 ) << 2 << 5 << 2 << -4 << endr
                              << -4 << -2 << -4 << -5 << endr
                              << -4 << -1 << 5 << 2 << endr;

        shifted = FFTHelper::fftShift(sig);
        for(int k = 0; k < int(sig.n_slices); k++){
            for(int i = 0; i < int(sig.n_rows); i++){
                for(int j = 0; j < int(sig.n_cols); j++){
                     CHECK_EQUAL(real(shifted(i,j,k)), real(centered(i,j,k)));
                }
             }
        }
        /*------------------------------------------------------*/

        sig.set_size( 5, 5, 2 );
        centered.set_size( 5, 5, 2 );
        sig.slice( 0 ) << -4 << 3 << -1 << 2 << -2 << endr
                        << 5 << 3 << 5 << -2 << -1 << endr
                        << 1 << 3 << -4 << 4 << -3 << endr
                        << -3 << 4 << 3 << 5 << 1 << endr
                        << 4 << 2 << 0 << -3 << -3 << endr;
        sig.slice( 1 ) << 3 << 2 << 3 << 3 << -5 << endr
                        << 5 << 4 << -1 << -4 << 3 << endr
                        << -4 << -5 << -4 << -5 << 3 << endr
                        << 0 << -4 << 2 << 0 << 0 << endr
                        << -5 << 5 << -4 << 0 << -4 << endr;

        centered.slice( 0 ) << 0 << 0 << 0 << -4 << 2 << endr
                              << 0 << -4 << -5 << 5 << -4 << endr
                              << 3 << -5 << 3 << 2 << 3 << endr
                              << -4 << 3 << 5 << 4 << -1 << endr
                              << -5 << 3 << -4 << -5 << -4 << endr;
        centered.slice( 1 ) << 5 << 1 << -3 << 4 << 3 << endr
                              << -3 << -3 << 4 << 2 << 0 << endr
                              << 2 << -2 << -4 << 3 << -1 << endr
                              << -2 << -1 << 5 << 3 << 5 << endr
                              << 4 << -3 << 1 << 3 << -4 << endr;

        shifted = FFTHelper::fftShift(sig);
        for(int k = 0; k < int(sig.n_slices); k++){
            for(int i = 0; i < int(sig.n_rows); i++){
                for(int j = 0; j < int(sig.n_cols); j++){
                     CHECK_EQUAL(real(shifted(i,j,k)), real(centered(i,j,k)));
                }
             }
        }
        /*------------------------------------------------------*/

        sig.set_size( 5, 3, 5 );
        centered.set_size( 5, 3, 5 );
        sig.slice( 0 ) << -2 << -1 << 0 << endr
                        << -1 << -4 << 5 << endr
                        << 4 << -3 << -5 << endr
                        << -1 << -2 << 0 << endr
                        << 5 << -4 << 0 << endr;
        sig.slice( 1 ) << 2 << 4 << 1 << endr
                        << -2 << -4 << 4 << endr
                        << 1 << 3 << -4 << endr
                        << -3 << -3 << -2 << endr
                        << -2 << 5 << -1 << endr;
        sig.slice( 2 ) << 4 << 5 << -2 << endr
                        << 3 << -3 << 5 << endr
                        << 3 << 1 << -2 << endr
                        << 0 << 3 << 2 << endr
                        << -1 << 1 << 2 << endr;
        sig.slice( 3 ) << -5 << -1 << 1 << endr
                        << -1 << 2 << -3 << endr
                        << -3 << -3 << -2 << endr
                        << 4 << -4 << 0 << endr
                        << -5 << 3 << 5 << endr;
        sig.slice( 4 ) << 3 << 3 << 2 << endr
                        << 5 << 4 << 3 << endr
                        << 3 << 0 << -3 << endr
                        << 3 << 2 << 0 << endr
                        << -1 << -3 << -4 << endr;

        centered.slice( 0 ) << 0 << 4 << -4 << endr
                              << 5 << -5 << 3 << endr
                              << 1 << -5 << -1 << endr
                              << -3 << -1 << 2 << endr
                              << -2 << -3 << -3 << endr;
        centered.slice( 1 ) << 0 << 3 << 2 << endr
                              << -4 << -1 << -3 << endr
                              << 2 << 3 << 3 << endr
                              << 3 << 5 << 4 << endr
                              << -3 << 3 << 0 << endr;
        centered.slice( 2 ) << 0 << -1 << -2 << endr
                              << 0 << 5 << -4 << endr
                              << 0 << -2 << -1 << endr
                              << 5 << -1 << -4 << endr
                              << -5 << 4 << -3 << endr;
        centered.slice( 3 ) << -2 << -3 << -3 << endr
                              << -1 << -2 << 5 << endr
                              << 1 << 2 << 4 << endr
                              << 4 << -2 << -4 << endr
                              << -4 << 1 << 3 << endr;
        centered.slice( 4 ) << 2 << 0 << 3 << endr
                              << 2 << -1 << 1 << endr
                              << -2 << 4 << 5 << endr
                              << 5 << 3 << -3 << endr
                              << -2 << 3 << 1 << endr;

        shifted = FFTHelper::fftShift(sig);
        for(int k = 0; k < int(sig.n_slices); k++){
            for(int i = 0; i < int(sig.n_rows); i++){
                for(int j = 0; j < int(sig.n_cols); j++){
                     CHECK_EQUAL(real(shifted(i,j,k)), real(centered(i,j,k)));
                }
             }
        }
        /*------------------------------------------------------*/

        sig.set_size( 2, 4, 2 );
        centered.set_size( 2, 4, 2 );
        sig.slice( 0 ) << 4 << 3 << -1 << -1 << endr
                        << 4 << 5 << 2 << -1 << endr;
        sig.slice( 1 ) << 2 << 4 << 2 << -4 << endr
                        << -1 << -2 << 5 << -4 << endr;

        centered.slice( 0 ) << 5 << -4 << -1 << -2 << endr
                              << 2 << -4 << 2 << 4 << endr;
        centered.slice( 1 ) << 2 << -1 << 4 << 5 << endr
                              << -1 << -1 << 4 << 3 << endr;

        shifted = FFTHelper::fftShift(sig);
        for(int k = 0; k < int(sig.n_slices); k++){
            for(int i = 0; i < int(sig.n_rows); i++){
                for(int j = 0; j < int(sig.n_cols); j++){
                     CHECK_EQUAL(real(shifted(i,j,k)), real(centered(i,j,k)));
                }
             }
        }
        /*------------------------------------------------------*/

        sig.set_size( 3, 5, 2 );
        centered.set_size( 3, 5, 2 );
        sig.slice( 0 ) << -1 << 4 << -4 << 1 << -4 << endr
                        << 3 << -4 << 1 << 5 << 0 << endr
                        << 5 << -5 << -3 << -3 << -2 << endr;
        sig.slice( 1 ) << 1 << 2 << 2 << -4 << 2 << endr
                        << 0 << 3 << 4 << 1 << -3 << endr
                        << 5 << 1 << 3 << -3 << -5 << endr;

        centered.slice( 0 ) << -3 << -5 << 5 << 1 << 3 << endr
                              << -4 << 2 << 1 << 2 << 2 << endr
                              << 1 << -3 << 0 << 3 << 4 << endr;
        centered.slice( 1 ) << -3 << -2 << 5 << -5 << -3 << endr
                              << 1 << -4 << -1 << 4 << -4 << endr
                              << 5 << 0 << 3 << -4 << 1 << endr;

        shifted = FFTHelper::fftShift(sig);
        for(int k = 0; k < int(sig.n_slices); k++){
            for(int i = 0; i < int(sig.n_rows); i++){
                for(int j = 0; j < int(sig.n_cols); j++){
                     CHECK_EQUAL(real(shifted(i,j,k)), real(centered(i,j,k)));
                }
             }
        }
        /*------------------------------------------------------*/

        sig.set_size( 5, 5, 4 );
        centered.set_size( 5, 5, 4 );
        sig.slice( 0 ) << -1 << -2 << 5 << -1 << 0 << endr
                        << -4 << -3 << -1 << 0 << 0 << endr
                        << -3 << -4 << 3 << 5 << -2 << endr
                        << -2 << 0 << 1 << -5 << -5 << endr
                        << -1 << 0 << -3 << -1 << 2 << endr;
        sig.slice( 1 ) << 5 << 3 << -4 << 4 << 2 << endr
                        << 5 << -4 << 2 << -3 << 4 << endr
                        << 3 << 2 << 3 << 4 << -1 << endr
                        << 1 << 3 << 5 << 0 << -2 << endr
                        << 1 << -5 << -2 << -4 << -5 << endr;
        sig.slice( 2 ) << -1 << -3 << 5 << 2 << 1 << endr
                        << 1 << 2 << -5 << -5 << 1 << endr
                        << -1 << 2 << 4 << 3 << 0 << endr
                        << -3 << -1 << 0 << 5 << 0 << endr
                        << -4 << -2 << 4 << 1 << 2 << endr;
        sig.slice( 3 ) << -2 << -3 << -3 << 3 << -1 << endr
                        << -1 << 1 << -1 << -5 << -3 << endr
                        << -1 << 2 << -3 << -1 << 2 << endr
                        << 4 << 4 << 5 << 3 << -5 << endr
                        << 5 << -5 << -2 << -2 << -1 << endr;

        centered.slice( 0 ) << 5 << 0 << -3 << -1 << 0 << endr
                              << 1 << 2 << -4 << -2 << 4 << endr
                              << 2 << 1 << -1 << -3 << 5 << endr
                              << -5 << 1 << 1 << 2 << -5 << endr
                              << 3 << 0 << -1 << 2 << 4 << endr;
        centered.slice( 1 ) << 3 << -5 << 4 << 4 << 5 << endr
                              << -2 << -1 << 5 << -5 << -2 << endr
                              << 3 << -1 << -2 << -3 << -3 << endr
                              << -5 << -3 << -1 << 1 << -1 << endr
                              << -1 << 2 << -1 << 2 << -3 << endr;
        centered.slice( 2 ) << -5 << -5 << -2 << 0 << 1 << endr
                              << -1 << 2 << -1 << 0 << -3 << endr
                              << -1 << 0 << -1 << -2 << 5 << endr
                              << 0 << 0 << -4 << -3 << -1 << endr
                              << 5 << -2 << -3 << -4 << 3 << endr;
        centered.slice( 3 ) << 0 << -2 << 1 << 3 << 5 << endr
                              << -4 << -5 << 1 << -5 << -2 << endr
                              << 4 << 2 << 5 << 3 << -4 << endr
                              << -3 << 4 << 5 << -4 << 2 << endr
                              << 4 << -1 << 3 << 2 << 3 << endr;

        shifted = FFTHelper::fftShift(sig);
        for(int k = 0; k < int(sig.n_slices); k++){
            for(int i = 0; i < int(sig.n_rows); i++){
                for(int j = 0; j < int(sig.n_cols); j++){
                     CHECK_EQUAL(real(shifted(i,j,k)), real(centered(i,j,k)));
                }
             }
        }
        /*------------------------------------------------------*/

        sig.set_size( 4, 4, 4 );
        centered.set_size( 4, 4, 4 );
        sig.slice( 0 ) << 5 << -1 << 5 << 3 << endr
                        << -3 << 5 << -2 << 0 << endr
                        << 1 << -2 << -5 << 1 << endr
                        << 1 << -4 << 1 << -4 << endr;
        sig.slice( 1 ) << -5 << -1 << 0 << 0 << endr
                        << -2 << -1 << 1 << -1 << endr
                        << 1 << 3 << 5 << 5 << endr
                        << 2 << 0 << -1 << -5 << endr;
        sig.slice( 2 ) << 2 << 1 << 4 << -3 << endr
                        << 2 << -5 << -4 << -2 << endr
                        << -4 << -2 << -1 << -1 << endr
                        << 2 << 2 << 3 << 5 << endr;
        sig.slice( 3 ) << -5 << -4 << -2 << -4 << endr
                        << 1 << -1 << 1 << 3 << endr
                        << -1 << 0 << 3 << 3 << endr
                        << -2 << -2 << -3 << 4 << endr;

        centered.slice( 0 ) << -1 << -1 << -4 << -2 << endr
                              << 3 << 5 << 2 << 2 << endr
                              << 4 << -3 << 2 << 1 << endr
                              << -4 << -2 << 2 << -5 << endr;
        centered.slice( 1 ) << 3 << 3 << -1 << 0 << endr
                              << -3 << 4 << -2 << -2 << endr
                              << -2 << -4 << -5 << -4 << endr
                              << 1 << 3 << 1 << -1 << endr;
        centered.slice( 2 ) << -5 << 1 << 1 << -2 << endr
                              << 1 << -4 << 1 << -4 << endr
                              << 5 << 3 << 5 << -1 << endr
                              << -2 << 0 << -3 << 5 << endr;
        centered.slice( 3 ) << 5 << 5 << 1 << 3 << endr
                              << -1 << -5 << 2 << 0 << endr
                              << 0 << 0 << -5 << -1 << endr
                              << 1 << -1 << -2 << -1 << endr;

        shifted = FFTHelper::fftShift(sig);
        for(int k = 0; k < int(sig.n_slices); k++){
            for(int i = 0; i < int(sig.n_rows); i++){
                for(int j = 0; j < int(sig.n_cols); j++){
                     CHECK_EQUAL(real(shifted(i,j,k)), real(centered(i,j,k)));
                }
             }
        }
        /*------------------------------------------------------*/

        sig.set_size( 5, 2, 2 );
        centered.set_size( 5, 2, 2 );
        sig.slice( 0 ) << 0 << 0 << endr
                        << 2 << -4 << endr
                        << -4 << 2 << endr
                        << -5 << 2 << endr
                        << 1 << -4 << endr;
        sig.slice( 1 ) << 3 << -1 << endr
                        << -2 << -1 << endr
                        << 2 << 3 << endr
                        << 3 << 5 << endr
                        << -2 << 4 << endr;

        centered.slice( 0 ) << 5 << 3 << endr
                              << 4 << -2 << endr
                              << -1 << 3 << endr
                              << -1 << -2 << endr
                              << 3 << 2 << endr;
        centered.slice( 1 ) << 2 << -5 << endr
                              << -4 << 1 << endr
                              << 0 << 0 << endr
                              << -4 << 2 << endr
                              << 2 << -4 << endr;

        shifted = FFTHelper::fftShift(sig);
        for(int k = 0; k < int(sig.n_slices); k++){
            for(int i = 0; i < int(sig.n_rows); i++){
                for(int j = 0; j < int(sig.n_cols); j++){
                     CHECK_EQUAL(real(shifted(i,j,k)), real(centered(i,j,k)));
                }
             }
        }
        /*------------------------------------------------------*/

        sig.set_size( 2, 2, 4 );
        centered.set_size( 2, 2, 4 );
        sig.slice( 0 ) << -1 << 5 << endr
                        << 3 << -2 << endr;
        sig.slice( 1 ) << -1 << 5 << endr
                        << -4 << -3 << endr;
        sig.slice( 2 ) << 2 << 2 << endr
                        << -2 << -1 << endr;
        sig.slice( 3 ) << 4 << -5 << endr
                        << 0 << 3 << endr;

        centered.slice( 0 ) << -1 << -2 << endr
                              << 2 << 2 << endr;
        centered.slice( 1 ) << 3 << 0 << endr
                              << -5 << 4 << endr;
        centered.slice( 2 ) << -2 << 3 << endr
                              << 5 << -1 << endr;
        centered.slice( 3 ) << -3 << -4 << endr
                              << 5 << -1 << endr;

        shifted = FFTHelper::fftShift(sig);
        for(int k = 0; k < int(sig.n_slices); k++){
            for(int i = 0; i < int(sig.n_rows); i++){
                for(int j = 0; j < int(sig.n_cols); j++){
                     CHECK_EQUAL(real(shifted(i,j,k)), real(centered(i,j,k)));
                }
             }
        }
        /*------------------------------------------------------*/

        sig.set_size( 5, 3, 4 );
        centered.set_size( 5, 3, 4 );
        sig.slice( 0 ) << 3 << -4 << 2 << endr
                        << 2 << -2 << -2 << endr
                        << -3 << 0 << 4 << endr
                        << 2 << 1 << -2 << endr
                        << -4 << -5 << 2 << endr;
        sig.slice( 1 ) << 3 << -2 << 5 << endr
                        << -3 << 2 << -1 << endr
                        << 0 << 3 << -1 << endr
                        << -5 << 4 << -2 << endr
                        << 1 << -3 << -5 << endr;
        sig.slice( 2 ) << -2 << 5 << 1 << endr
                        << 3 << 0 << 4 << endr
                        << 2 << -2 << -4 << endr
                        << 0 << 4 << -1 << endr
                        << -5 << 2 << 5 << endr;
        sig.slice( 3 ) << 5 << 3 << 3 << endr
                        << 1 << -5 << -2 << endr
                        << 0 << 1 << -4 << endr
                        << -3 << 4 << -5 << endr
                        << 4 << 3 << 3 << endr;

        centered.slice( 0 ) << -1 << 0 << 4 << endr
                              << 5 << -5 << 2 << endr
                              << 1 << -2 << 5 << endr
                              << 4 << 3 << 0 << endr
                              << -4 << 2 << -2 << endr;
        centered.slice( 1 ) << -5 << -3 << 4 << endr
                              << 3 << 4 << 3 << endr
                              << 3 << 5 << 3 << endr
                              << -2 << 1 << -5 << endr
                              << -4 << 0 << 1 << endr;
        centered.slice( 2 ) << -2 << 2 << 1 << endr
                              << 2 << -4 << -5 << endr
                              << 2 << 3 << -4 << endr
                              << -2 << 2 << -2 << endr
                              << 4 << -3 << 0 << endr;
        centered.slice( 3 ) << -2 << -5 << 4 << endr
                              << -5 << 1 << -3 << endr
                              << 5 << 3 << -2 << endr
                              << -1 << -3 << 2 << endr
                              << -1 << 0 << 3 << endr;

        shifted = FFTHelper::fftShift(sig);
        for(int k = 0; k < int(sig.n_slices); k++){
            for(int i = 0; i < int(sig.n_rows); i++){
                for(int j = 0; j < int(sig.n_cols); j++){
                     CHECK_EQUAL(real(shifted(i,j,k)), real(centered(i,j,k)));
                }
             }
        }
        /*------------------------------------------------------*/

        sig.set_size( 4, 3, 4 );
        centered.set_size( 4, 3, 4 );
        sig.slice( 0 ) << 1 << -4 << 1 << endr
                        << -1 << 5 << 0 << endr
                        << 3 << -3 << -2 << endr
                        << -2 << -4 << 4 << endr;
        sig.slice( 1 ) << -3 << -4 << 4 << endr
                        << 0 << -3 << 0 << endr
                        << 3 << 0 << -4 << endr
                        << 0 << 3 << 5 << endr;
        sig.slice( 2 ) << -5 << 1 << -4 << endr
                        << -3 << -4 << -4 << endr
                        << -5 << -3 << 1 << endr
                        << 5 << -5 << -4 << endr;
        sig.slice( 3 ) << -5 << -4 << -1 << endr
                        << 3 << 3 << -3 << endr
                        << -2 << 5 << -3 << endr
                        << -2 << -4 << 3 << endr;

        centered.slice( 0 ) << 1 << -5 << -3 << endr
                              << -4 << 5 << -5 << endr
                              << -4 << -5 << 1 << endr
                              << -4 << -3 << -4 << endr;
        centered.slice( 1 ) << -3 << -2 << 5 << endr
                              << 3 << -2 << -4 << endr
                              << -1 << -5 << -4 << endr
                              << -3 << 3 << 3 << endr;
        centered.slice( 2 ) << -2 << 3 << -3 << endr
                              << 4 << -2 << -4 << endr
                              << 1 << 1 << -4 << endr
                              << 0 << -1 << 5 << endr;
        centered.slice( 3 ) << -4 << 3 << 0 << endr
                              << 5 << 0 << 3 << endr
                              << 4 << -3 << -4 << endr
                              << 0 << 0 << -3 << endr;

        shifted = FFTHelper::fftShift(sig);
        for(int k = 0; k < int(sig.n_slices); k++){
            for(int i = 0; i < int(sig.n_rows); i++){
                for(int j = 0; j < int(sig.n_cols); j++){
                     CHECK_EQUAL(real(shifted(i,j,k)), real(centered(i,j,k)));
                }
             }
        }
        /*------------------------------------------------------*/

        sig.set_size( 2, 5, 3 );
        centered.set_size( 2, 5, 3 );
        sig.slice( 0 ) << 4 << -5 << 1 << 0 << -4 << endr
                        << -5 << -5 << 1 << 2 << -1 << endr;
        sig.slice( 1 ) << 3 << 0 << 1 << 4 << 0 << endr
                        << -2 << 5 << 3 << -1 << 5 << endr;
        sig.slice( 2 ) << -2 << -4 << -3 << -5 << -3 << endr
                        << -2 << 1 << 0 << -3 << -5 << endr;

        centered.slice( 0 ) << -3 << -5 << -2 << 1 << 0 << endr
                              << -5 << -3 << -2 << -4 << -3 << endr;
        centered.slice( 1 ) << 2 << -1 << -5 << -5 << 1 << endr
                              << 0 << -4 << 4 << -5 << 1 << endr;
        centered.slice( 2 ) << -1 << 5 << -2 << 5 << 3 << endr
                              << 4 << 0 << 3 << 0 << 1 << endr;

        shifted = FFTHelper::fftShift(sig);
        for(int k = 0; k < int(sig.n_slices); k++){
            for(int i = 0; i < int(sig.n_rows); i++){
                for(int j = 0; j < int(sig.n_cols); j++){
                     CHECK_EQUAL(real(shifted(i,j,k)), real(centered(i,j,k)));
                }
             }
        }
        /*------------------------------------------------------*/

        sig.set_size( 2, 3, 4 );
        centered.set_size( 2, 3, 4 );
        sig.slice( 0 ) << 3 << 0 << -1 << endr
                        << -5 << -2 << -5 << endr;
        sig.slice( 1 ) << 1 << 3 << 0 << endr
                        << -4 << 0 << 2 << endr;
        sig.slice( 2 ) << 0 << -2 << 0 << endr
                        << 5 << 2 << 3 << endr;
        sig.slice( 3 ) << 5 << 5 << -4 << endr
                        << 4 << 3 << 5 << endr;

        centered.slice( 0 ) << 3 << 5 << 2 << endr
                              << 0 << 0 << -2 << endr;
        centered.slice( 1 ) << 5 << 4 << 3 << endr
                              << -4 << 5 << 5 << endr;
        centered.slice( 2 ) << -5 << -5 << -2 << endr
                              << -1 << 3 << 0 << endr;
        centered.slice( 3 ) << 2 << -4 << 0 << endr
                              << 0 << 1 << 3 << endr;

        shifted = FFTHelper::fftShift(sig);
        for(int k = 0; k < int(sig.n_slices); k++){
            for(int i = 0; i < int(sig.n_rows); i++){
                for(int j = 0; j < int(sig.n_cols); j++){
                     CHECK_EQUAL(real(shifted(i,j,k)), real(centered(i,j,k)));
                }
             }
        }
        /*------------------------------------------------------*/

        sig.set_size( 4, 4, 5 );
        centered.set_size( 4, 4, 5 );
        sig.slice( 0 ) << 3 << 3 << -3 << 1 << endr
                        << -1 << -2 << 1 << 4 << endr
                        << -2 << -2 << 5 << -1 << endr
                        << 4 << 2 << 3 << -2 << endr;
        sig.slice( 1 ) << 0 << 2 << 4 << 1 << endr
                        << 3 << 4 << 4 << -3 << endr
                        << 0 << -1 << 5 << -2 << endr
                        << 3 << -4 << -2 << 4 << endr;
        sig.slice( 2 ) << 0 << -4 << 2 << -5 << endr
                        << -4 << 5 << 2 << -3 << endr
                        << -4 << -3 << -4 << -1 << endr
                        << 1 << 4 << -3 << 1 << endr;
        sig.slice( 3 ) << -4 << -3 << 1 << 3 << endr
                        << 1 << 4 << -1 << 4 << endr
                        << -1 << 4 << -2 << 0 << endr
                        << -3 << -4 << -3 << 2 << endr;
        sig.slice( 4 ) << 0 << -3 << 2 << 1 << endr
                        << -1 << 2 << -3 << -1 << endr
                        << -4 << 2 << -5 << 3 << endr
                        << 1 << -4 << 3 << -5 << endr;

        centered.slice( 0 ) << -2 << 0 << -1 << 4 << endr
                              << -3 << 2 << -3 << -4 << endr
                              << 1 << 3 << -4 << -3 << endr
                              << -1 << 4 << 1 << 4 << endr;
        centered.slice( 1 ) << -5 << 3 << -4 << 2 << endr
                              << 3 << -5 << 1 << -4 << endr
                              << 2 << 1 << 0 << -3 << endr
                              << -3 << -1 << -1 << 2 << endr;
        centered.slice( 2 ) << 5 << -1 << -2 << -2 << endr
                              << 3 << -2 << 4 << 2 << endr
                              << -3 << 1 << 3 << 3 << endr
                              << 1 << 4 << -1 << -2 << endr;
        centered.slice( 3 ) << 5 << -2 << 0 << -1 << endr
                              << -2 << 4 << 3 << -4 << endr
                              << 4 << 1 << 0 << 2 << endr
                              << 4 << -3 << 3 << 4 << endr;
        centered.slice( 4 ) << -4 << -1 << -4 << -3 << endr
                              << -3 << 1 << 1 << 4 << endr
                              << 2 << -5 << 0 << -4 << endr
                              << 2 << -3 << -4 << 5 << endr;

        shifted = FFTHelper::fftShift(sig);
        for(int k = 0; k < int(sig.n_slices); k++){
            for(int i = 0; i < int(sig.n_rows); i++){
                for(int j = 0; j < int(sig.n_cols); j++){
                     CHECK_EQUAL(real(shifted(i,j,k)), real(centered(i,j,k)));
                }
             }
        }
        /*------------------------------------------------------*/

        sig.set_size( 4, 4, 2 );
        centered.set_size( 4, 4, 2 );
        sig.slice( 0 ) << -5 << -5 << 5 << -5 << endr
                        << 3 << -2 << -3 << -1 << endr
                        << 2 << -5 << -5 << 5 << endr
                        << -2 << -1 << 5 << 4 << endr;
        sig.slice( 1 ) << -5 << 3 << 1 << -4 << endr
                        << 0 << -4 << 3 << -4 << endr
                        << -5 << -4 << 5 << -3 << endr
                        << 0 << -1 << -2 << -4 << endr;

        centered.slice( 0 ) << 5 << -3 << -5 << -4 << endr
                              << -2 << -4 << 0 << -1 << endr
                              << 1 << -4 << -5 << 3 << endr
                              << 3 << -4 << 0 << -4 << endr;
        centered.slice( 1 ) << -5 << 5 << 2 << -5 << endr
                              << 5 << 4 << -2 << -1 << endr
                              << 5 << -5 << -5 << -5 << endr
                              << -3 << -1 << 3 << -2 << endr;

        shifted = FFTHelper::fftShift(sig);
        for(int k = 0; k < int(sig.n_slices); k++){
            for(int i = 0; i < int(sig.n_rows); i++){
                for(int j = 0; j < int(sig.n_cols); j++){
                     CHECK_EQUAL(real(shifted(i,j,k)), real(centered(i,j,k)));
                }
             }
        }
        /*------------------------------------------------------*/

        sig.set_size( 5, 4, 5 );
        centered.set_size( 5, 4, 5 );
        sig.slice( 0 ) << -3 << -1 << 4 << 4 << endr
                        << 4 << 5 << -1 << -3 << endr
                        << 0 << -1 << -1 << 3 << endr
                        << 3 << -2 << 3 << -1 << endr
                        << 5 << -3 << 5 << 2 << endr;
        sig.slice( 1 ) << -3 << -5 << -3 << -1 << endr
                        << 3 << 0 << 5 << 3 << endr
                        << -5 << 1 << -1 << -5 << endr
                        << 5 << 4 << 5 << 1 << endr
                        << 2 << -3 << -2 << 5 << endr;
        sig.slice( 2 ) << -4 << 2 << 3 << 5 << endr
                        << 2 << -5 << 0 << 4 << endr
                        << 5 << -5 << -2 << -1 << endr
                        << 2 << 2 << 0 << 0 << endr
                        << 5 << -3 << 0 << 2 << endr;
        sig.slice( 3 ) << -2 << 4 << 1 << 2 << endr
                        << -3 << -3 << -3 << -1 << endr
                        << 5 << 0 << -1 << -5 << endr
                        << 5 << 5 << -4 << -1 << endr
                        << 2 << 0 << -2 << -4 << endr;
        sig.slice( 4 ) << -5 << -4 << 2 << 0 << endr
                        << -5 << 3 << 0 << 5 << endr
                        << 0 << -4 << -4 << 5 << endr
                        << -4 << 2 << -2 << 2 << endr
                        << -4 << 1 << -4 << -5 << endr;

        centered.slice( 0 ) << -4 << -1 << 5 << 5 << endr
                              << -2 << -4 << 2 << 0 << endr
                              << 1 << 2 << -2 << 4 << endr
                              << -3 << -1 << -3 << -3 << endr
                              << -1 << -5 << 5 << 0 << endr;
        centered.slice( 1 ) << -2 << 2 << -4 << 2 << endr
                              << -4 << -5 << -4 << 1 << endr
                              << 2 << 0 << -5 << -4 << endr
                              << 0 << 5 << -5 << 3 << endr
                              << -4 << 5 << 0 << -4 << endr;
        centered.slice( 2 ) << 3 << -1 << 3 << -2 << endr
                              << 5 << 2 << 5 << -3 << endr
                              << 4 << 4 << -3 << -1 << endr
                              << -1 << -3 << 4 << 5 << endr
                              << -1 << 3 << 0 << -1 << endr;
        centered.slice( 3 ) << 5 << 1 << 5 << 4 << endr
                              << -2 << 5 << 2 << -3 << endr
                              << -3 << -1 << -3 << -5 << endr
                              << 5 << 3 << 3 << 0 << endr
                              << -1 << -5 << -5 << 1 << endr;
        centered.slice( 4 ) << 0 << 0 << 2 << 2 << endr
                              << 0 << 2 << 5 << -3 << endr
                              << 3 << 5 << -4 << 2 << endr
                              << 0 << 4 << 2 << -5 << endr
                              << -2 << -1 << 5 << -5 << endr;

        shifted = FFTHelper::fftShift(sig);
        for(int k = 0; k < int(sig.n_slices); k++){
            for(int i = 0; i < int(sig.n_rows); i++){
                for(int j = 0; j < int(sig.n_cols); j++){
                     CHECK_EQUAL(real(shifted(i,j,k)), real(centered(i,j,k)));
                }
             }
        }
        /*------------------------------------------------------*/

        sig.set_size( 5, 4, 5 );
        centered.set_size( 5, 4, 5 );
        sig.slice( 0 ) << -2 << 4 << 2 << -4 << endr
                        << 3 << 2 << 3 << 4 << endr
                        << 4 << -3 << -2 << -5 << endr
                        << 4 << 0 << 2 << -1 << endr
                        << -2 << -5 << 1 << 4 << endr;
        sig.slice( 1 ) << -1 << 0 << 2 << 4 << endr
                        << -4 << 5 << 1 << -4 << endr
                        << 2 << 3 << 4 << -5 << endr
                        << 3 << -4 << -3 << -2 << endr
                        << 1 << 5 << -5 << 3 << endr;
        sig.slice( 2 ) << -3 << -4 << 2 << 5 << endr
                        << 4 << 0 << -2 << 2 << endr
                        << -3 << -1 << -5 << -2 << endr
                        << 1 << -2 << 3 << 2 << endr
                        << 4 << -4 << -2 << 1 << endr;
        sig.slice( 3 ) << 4 << -3 << -1 << -3 << endr
                        << -1 << -3 << 0 << -5 << endr
                        << 4 << 1 << 3 << 2 << endr
                        << -3 << 4 << 2 << 4 << endr
                        << 4 << 2 << -4 << 1 << endr;
        sig.slice( 4 ) << -1 << -2 << 4 << 2 << endr
                        << -2 << 3 << 0 << -5 << endr
                        << -4 << -3 << -3 << 2 << endr
                        << 0 << -5 << -2 << 1 << endr
                        << 1 << 4 << -5 << 2 << endr;

        centered.slice( 0 ) << 2 << 4 << -3 << 4 << endr
                              << -4 << 1 << 4 << 2 << endr
                              << -1 << -3 << 4 << -3 << endr
                              << 0 << -5 << -1 << -3 << endr
                              << 3 << 2 << 4 << 1 << endr;
        centered.slice( 1 ) << -2 << 1 << 0 << -5 << endr
                              << -5 << 2 << 1 << 4 << endr
                              << 4 << 2 << -1 << -2 << endr
                              << 0 << -5 << -2 << 3 << endr
                              << -3 << 2 << -4 << -3 << endr;
        centered.slice( 2 ) << 2 << -1 << 4 << 0 << endr
                              << 1 << 4 << -2 << -5 << endr
                              << 2 << -4 << -2 << 4 << endr
                              << 3 << 4 << 3 << 2 << endr
                              << -2 << -5 << 4 << -3 << endr;
        centered.slice( 3 ) << -3 << -2 << 3 << -4 << endr
                              << -5 << 3 << 1 << 5 << endr
                              << 2 << 4 << -1 << 0 << endr
                              << 1 << -4 << -4 << 5 << endr
                              << 4 << -5 << 2 << 3 << endr;
        centered.slice( 4 ) << 3 << 2 << 1 << -2 << endr
                              << -2 << 1 << 4 << -4 << endr
                              << 2 << 5 << -3 << -4 << endr
                              << -2 << 2 << 4 << 0 << endr
                              << -5 << -2 << -3 << -1 << endr;

        shifted = FFTHelper::fftShift(sig);
        for(int k = 0; k < int(sig.n_slices); k++){
            for(int i = 0; i < int(sig.n_rows); i++){
                for(int j = 0; j < int(sig.n_cols); j++){
                     CHECK_EQUAL(real(shifted(i,j,k)), real(centered(i,j,k)));
                }
             }
        }
        /*------------------------------------------------------*/

        sig.set_size( 5, 3, 2 );
        centered.set_size( 5, 3, 2 );
        sig.slice( 0 ) << -3 << -2 << -4 << endr
                        << -3 << -1 << -5 << endr
                        << -4 << -2 << -4 << endr
                        << 0 << -2 << 0 << endr
                        << 0 << 2 << 5 << endr;
        sig.slice( 1 ) << -2 << 4 << 1 << endr
                        << 0 << -3 << 0 << endr
                        << 1 << 1 << -2 << endr
                        << -3 << 5 << 2 << endr
                        << -4 << -2 << 3 << endr;

        centered.slice( 0 ) << 2 << -3 << 5 << endr
                              << 3 << -4 << -2 << endr
                              << 1 << -2 << 4 << endr
                              << 0 << 0 << -3 << endr
                              << -2 << 1 << 1 << endr;
        centered.slice( 1 ) << 0 << 0 << -2 << endr
                              << 5 << 0 << 2 << endr
                              << -4 << -3 << -2 << endr
                              << -5 << -3 << -1 << endr
                              << -4 << -4 << -2 << endr;

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
