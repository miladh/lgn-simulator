/**********************************************************************
 *  Test: fftShift function (2d):
 *        shift the zero-frequency component to the center of the spectrum.
 *
 *  Analytic source: Pyhton's numpy.fft.fftshift
 *
 * ********************************************************************/

#include <unittest++/UnitTest++.h>
#include <lgnSimulator.h>

using namespace lgnSimulator;

SUITE(fftHelper){


    TEST(fftshift_2D){
        cx_mat sig, centered, shifted;


        sig << 0 << -2 << -5 << endr
           << -1 << -4 << -5 << endr
           << -5 << -1 << 5 << endr
           << 1 << -4 << -1 << endr;

        centered << 5 << -5 << -1 << endr
                << -1 << 1 << -4 << endr
                << -5 << 0 << -2 << endr
                << -5 << -1 << -4 << endr;
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_rows); i++){
           for(int j = 0; j < int(sig.n_cols); j++){
            CHECK_EQUAL(real(shifted(i,j)), real(centered(i,j)));
            }
        }
        /*------------------------------------------------------*/

        sig << 3 << -3 << endr
           << 3 << -1 << endr;

        centered << -1 << 3 << endr
                << -3 << 3 << endr;
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_rows); i++){
           for(int j = 0; j < int(sig.n_cols); j++){
            CHECK_EQUAL(real(shifted(i,j)), real(centered(i,j)));
            }
        }
        /*------------------------------------------------------*/

        sig << 2 << -5 << 4 << -1 << -2 << endr
           << -3 << 2 << -2 << 0 << 0 << endr
           << 2 << 2 << 3 << 1 << 2 << endr
           << 1 << -1 << -4 << 5 << -1 << endr
           << 1 << -3 << 4 << 2 << 1 << endr;

        centered << 5 << -1 << 1 << -1 << -4 << endr
                << 2 << 1 << 1 << -3 << 4 << endr
                << -1 << -2 << 2 << -5 << 4 << endr
                << 0 << 0 << -3 << 2 << -2 << endr
                << 1 << 2 << 2 << 2 << 3 << endr;
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_rows); i++){
           for(int j = 0; j < int(sig.n_cols); j++){
            CHECK_EQUAL(real(shifted(i,j)), real(centered(i,j)));
            }
        }
        /*------------------------------------------------------*/

        sig << 1 << -4 << 1 << -3 << 1 << endr
           << 5 << 3 << -2 << 5 << 1 << endr;

        centered << 5 << 1 << 5 << 3 << -2 << endr
                << -3 << 1 << 1 << -4 << 1 << endr;
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_rows); i++){
           for(int j = 0; j < int(sig.n_cols); j++){
            CHECK_EQUAL(real(shifted(i,j)), real(centered(i,j)));
            }
        }
        /*------------------------------------------------------*/

        sig << -2 << 3 << -4 << 2 << -4 << endr
           << 4 << -3 << 2 << 2 << 4 << endr
           << 0 << 4 << -4 << -5 << -2 << endr
           << -4 << 3 << 1 << -1 << -2 << endr
           << -3 << 0 << 2 << 0 << -1 << endr;

        centered << -1 << -2 << -4 << 3 << 1 << endr
                << 0 << -1 << -3 << 0 << 2 << endr
                << 2 << -4 << -2 << 3 << -4 << endr
                << 2 << 4 << 4 << -3 << 2 << endr
                << -5 << -2 << 0 << 4 << -4 << endr;
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_rows); i++){
           for(int j = 0; j < int(sig.n_cols); j++){
            CHECK_EQUAL(real(shifted(i,j)), real(centered(i,j)));
            }
        }
        /*------------------------------------------------------*/

        sig << -2 << -1 << endr
           << 1 << 3 << endr
           << 4 << -1 << endr
           << -3 << -2 << endr
           << -3 << -1 << endr;

        centered << -2 << -3 << endr
                << -1 << -3 << endr
                << -1 << -2 << endr
                << 3 << 1 << endr
                << -1 << 4 << endr;
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_rows); i++){
           for(int j = 0; j < int(sig.n_cols); j++){
            CHECK_EQUAL(real(shifted(i,j)), real(centered(i,j)));
            }
        }
        /*------------------------------------------------------*/

        sig << -2 << 1 << -4 << endr
           << 5 << 5 << 4 << endr;

        centered << 4 << 5 << 5 << endr
                << -4 << -2 << 1 << endr;
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_rows); i++){
           for(int j = 0; j < int(sig.n_cols); j++){
            CHECK_EQUAL(real(shifted(i,j)), real(centered(i,j)));
            }
        }
        /*------------------------------------------------------*/

        sig << 5 << 5 << 1 << -3 << 5 << endr
           << 0 << 5 << -1 << 2 << 5 << endr
           << -2 << 1 << 4 << -4 << 0 << endr;

        centered << -4 << 0 << -2 << 1 << 4 << endr
                << -3 << 5 << 5 << 5 << 1 << endr
                << 2 << 5 << 0 << 5 << -1 << endr;
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_rows); i++){
           for(int j = 0; j < int(sig.n_cols); j++){
            CHECK_EQUAL(real(shifted(i,j)), real(centered(i,j)));
            }
        }
        /*------------------------------------------------------*/

        sig << 4 << 4 << endr
           << -2 << 3 << endr;

        centered << 3 << -2 << endr
                << 4 << 4 << endr;
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_rows); i++){
           for(int j = 0; j < int(sig.n_cols); j++){
            CHECK_EQUAL(real(shifted(i,j)), real(centered(i,j)));
            }
        }
        /*------------------------------------------------------*/

        sig << 0 << 1 << 4 << 1 << endr
           << 3 << 4 << -5 << -4 << endr
           << -3 << -5 << -5 << -5 << endr
           << 4 << -5 << 5 << -5 << endr
           << -1 << 0 << -4 << 4 << endr;

        centered << 5 << -5 << 4 << -5 << endr
                << -4 << 4 << -1 << 0 << endr
                << 4 << 1 << 0 << 1 << endr
                << -5 << -4 << 3 << 4 << endr
                << -5 << -5 << -3 << -5 << endr;
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_rows); i++){
           for(int j = 0; j < int(sig.n_cols); j++){
            CHECK_EQUAL(real(shifted(i,j)), real(centered(i,j)));
            }
        }
        /*------------------------------------------------------*/

        sig << -3 << -4 << -1 << 0 << endr
           << -1 << -4 << 0 << -1 << endr;

        centered << 0 << -1 << -1 << -4 << endr
                << -1 << 0 << -3 << -4 << endr;
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_rows); i++){
           for(int j = 0; j < int(sig.n_cols); j++){
            CHECK_EQUAL(real(shifted(i,j)), real(centered(i,j)));
            }
        }
        /*------------------------------------------------------*/

        sig << -2 << 2 << 5 << endr
           << -2 << 4 << 5 << endr
           << 5 << 4 << -5 << endr;

        centered << -5 << 5 << 4 << endr
                << 5 << -2 << 2 << endr
                << 5 << -2 << 4 << endr;
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_rows); i++){
           for(int j = 0; j < int(sig.n_cols); j++){
            CHECK_EQUAL(real(shifted(i,j)), real(centered(i,j)));
            }
        }
        /*------------------------------------------------------*/

        sig << 0 << 0 << 0 << 0 << endr
           << 1 << 1 << 4 << -4 << endr
           << 3 << -1 << 5 << -3 << endr
           << 0 << 5 << -5 << 3 << endr
           << 2 << 1 << 2 << 5 << endr;

        centered << -5 << 3 << 0 << 5 << endr
                << 2 << 5 << 2 << 1 << endr
                << 0 << 0 << 0 << 0 << endr
                << 4 << -4 << 1 << 1 << endr
                << 5 << -3 << 3 << -1 << endr;
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_rows); i++){
           for(int j = 0; j < int(sig.n_cols); j++){
            CHECK_EQUAL(real(shifted(i,j)), real(centered(i,j)));
            }
        }
        /*------------------------------------------------------*/

        sig << 0 << 4 << 1 << 2 << 5 << endr
           << 0 << -3 << -3 << 4 << 4 << endr;

        centered << 4 << 4 << 0 << -3 << -3 << endr
                << 2 << 5 << 0 << 4 << 1 << endr;
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_rows); i++){
           for(int j = 0; j < int(sig.n_cols); j++){
            CHECK_EQUAL(real(shifted(i,j)), real(centered(i,j)));
            }
        }
        /*------------------------------------------------------*/

        sig << 2 << 0 << -4 << 4 << 2 << endr
           << -5 << -2 << 1 << 1 << -1 << endr
           << 4 << 3 << 0 << -4 << 4 << endr
           << 0 << -2 << 4 << -3 << -2 << endr
           << 3 << 0 << -4 << 0 << -1 << endr;

        centered << -3 << -2 << 0 << -2 << 4 << endr
                << 0 << -1 << 3 << 0 << -4 << endr
                << 4 << 2 << 2 << 0 << -4 << endr
                << 1 << -1 << -5 << -2 << 1 << endr
                << -4 << 4 << 4 << 3 << 0 << endr;
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_rows); i++){
           for(int j = 0; j < int(sig.n_cols); j++){
            CHECK_EQUAL(real(shifted(i,j)), real(centered(i,j)));
            }
        }
        /*------------------------------------------------------*/

        sig << -4 << 3 << 2 << endr
           << -1 << -1 << 1 << endr;

        centered << 1 << -1 << -1 << endr
                << 2 << -4 << 3 << endr;
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_rows); i++){
           for(int j = 0; j < int(sig.n_cols); j++){
            CHECK_EQUAL(real(shifted(i,j)), real(centered(i,j)));
            }
        }
        /*------------------------------------------------------*/

        sig << -1 << 3 << 4 << 4 << endr
           << -2 << 5 << 1 << -1 << endr
           << -3 << -3 << -4 << 1 << endr;

        centered << -4 << 1 << -3 << -3 << endr
                << 4 << 4 << -1 << 3 << endr
                << 1 << -1 << -2 << 5 << endr;
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_rows); i++){
           for(int j = 0; j < int(sig.n_cols); j++){
            CHECK_EQUAL(real(shifted(i,j)), real(centered(i,j)));
            }
        }
        /*------------------------------------------------------*/

        sig << -2 << 3 << 4 << endr
           << -1 << 1 << 1 << endr
           << 1 << 3 << 3 << endr
           << -1 << -5 << -5 << endr;

        centered << 3 << 1 << 3 << endr
                << -5 << -1 << -5 << endr
                << 4 << -2 << 3 << endr
                << 1 << -1 << 1 << endr;
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_rows); i++){
           for(int j = 0; j < int(sig.n_cols); j++){
            CHECK_EQUAL(real(shifted(i,j)), real(centered(i,j)));
            }
        }
        /*------------------------------------------------------*/

        sig << 5 << 2 << -3 << endr
           << 3 << 1 << 5 << endr;

        centered << 5 << 3 << 1 << endr
                << -3 << 5 << 2 << endr;
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_rows); i++){
           for(int j = 0; j < int(sig.n_cols); j++){
            CHECK_EQUAL(real(shifted(i,j)), real(centered(i,j)));
            }
        }
        /*------------------------------------------------------*/

        sig << -4 << -5 << 2 << 1 << endr
           << -2 << -3 << 1 << 4 << endr;

        centered << 1 << 4 << -2 << -3 << endr
                << 2 << 1 << -4 << -5 << endr;
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_rows); i++){
           for(int j = 0; j < int(sig.n_cols); j++){
            CHECK_EQUAL(real(shifted(i,j)), real(centered(i,j)));
            }
        }
        /*------------------------------------------------------*/

        sig << -4 << -3 << 3 << endr
           << 0 << -4 << 2 << endr
           << -2 << -1 << 4 << endr;

        centered << 4 << -2 << -1 << endr
                << 3 << -4 << -3 << endr
                << 2 << 0 << -4 << endr;
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_rows); i++){
           for(int j = 0; j < int(sig.n_cols); j++){
            CHECK_EQUAL(real(shifted(i,j)), real(centered(i,j)));
            }
        }
        /*------------------------------------------------------*/

        sig << 5 << 5 << -1 << 4 << -2 << endr
           << -2 << -3 << -2 << 1 << -2 << endr;

        centered << 1 << -2 << -2 << -3 << -2 << endr
                << 4 << -2 << 5 << 5 << -1 << endr;
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_rows); i++){
           for(int j = 0; j < int(sig.n_cols); j++){
            CHECK_EQUAL(real(shifted(i,j)), real(centered(i,j)));
            }
        }
        /*------------------------------------------------------*/

        sig << -2 << 4 << -3 << -4 << endr
           << 3 << 3 << -4 << 3 << endr
           << -1 << 0 << 3 << -4 << endr
           << 4 << -4 << -2 << -1 << endr
           << 5 << -3 << -5 << 0 << endr;

        centered << -2 << -1 << 4 << -4 << endr
                << -5 << 0 << 5 << -3 << endr
                << -3 << -4 << -2 << 4 << endr
                << -4 << 3 << 3 << 3 << endr
                << 3 << -4 << -1 << 0 << endr;
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_rows); i++){
           for(int j = 0; j < int(sig.n_cols); j++){
            CHECK_EQUAL(real(shifted(i,j)), real(centered(i,j)));
            }
        }
        /*------------------------------------------------------*/

        sig << -1 << -5 << 4 << endr
           << -4 << 2 << 4 << endr
           << -1 << -4 << -4 << endr;

        centered << -4 << -1 << -4 << endr
                << 4 << -1 << -5 << endr
                << 4 << -4 << 2 << endr;
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_rows); i++){
           for(int j = 0; j < int(sig.n_cols); j++){
            CHECK_EQUAL(real(shifted(i,j)), real(centered(i,j)));
            }
        }
        /*------------------------------------------------------*/

        sig << -1 << -2 << -3 << -1 << -5 << endr
           << -1 << 5 << 5 << -4 << 5 << endr
           << 4 << -2 << 1 << -2 << 0 << endr
           << -5 << 5 << -4 << -2 << 1 << endr
           << -4 << -2 << -5 << -4 << -2 << endr;

        centered << -2 << 1 << -5 << 5 << -4 << endr
                << -4 << -2 << -4 << -2 << -5 << endr
                << -1 << -5 << -1 << -2 << -3 << endr
                << -4 << 5 << -1 << 5 << 5 << endr
                << -2 << 0 << 4 << -2 << 1 << endr;
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_rows); i++){
           for(int j = 0; j < int(sig.n_cols); j++){
            CHECK_EQUAL(real(shifted(i,j)), real(centered(i,j)));
            }
        }
        /*------------------------------------------------------*/

        sig << 2 << -2 << 1 << endr
           << -2 << -4 << -1 << endr
           << -5 << -4 << -1 << endr;

        centered << -1 << -5 << -4 << endr
                << 1 << 2 << -2 << endr
                << -1 << -2 << -4 << endr;
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_rows); i++){
           for(int j = 0; j < int(sig.n_cols); j++){
            CHECK_EQUAL(real(shifted(i,j)), real(centered(i,j)));
            }
        }
        /*------------------------------------------------------*/

        sig << 5 << 4 << endr
           << -4 << -5 << endr;

        centered << -5 << -4 << endr
                << 4 << 5 << endr;
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_rows); i++){
           for(int j = 0; j < int(sig.n_cols); j++){
            CHECK_EQUAL(real(shifted(i,j)), real(centered(i,j)));
            }
        }
        /*------------------------------------------------------*/

        sig << 1 << 4 << -1 << endr
           << -3 << -5 << -1 << endr
           << -1 << 4 << 3 << endr;

        centered << 3 << -1 << 4 << endr
                << -1 << 1 << 4 << endr
                << -1 << -3 << -5 << endr;
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_rows); i++){
           for(int j = 0; j < int(sig.n_cols); j++){
            CHECK_EQUAL(real(shifted(i,j)), real(centered(i,j)));
            }
        }
        /*------------------------------------------------------*/

        sig << -4 << -5 << endr
           << 5 << 1 << endr
           << -2 << 4 << endr
           << -1 << -2 << endr;

        centered << 4 << -2 << endr
                << -2 << -1 << endr
                << -5 << -4 << endr
                << 1 << 5 << endr;
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_rows); i++){
           for(int j = 0; j < int(sig.n_cols); j++){
            CHECK_EQUAL(real(shifted(i,j)), real(centered(i,j)));
            }
        }
        /*------------------------------------------------------*/

        sig << 3 << 5 << -1 << 3 << -2 << endr
           << 2 << 5 << 5 << -1 << 4 << endr
           << 1 << 1 << -3 << -2 << 2 << endr
           << -3 << -1 << -3 << -5 << 1 << endr
           << -4 << -4 << 5 << -2 << 2 << endr;

        centered << -5 << 1 << -3 << -1 << -3 << endr
                << -2 << 2 << -4 << -4 << 5 << endr
                << 3 << -2 << 3 << 5 << -1 << endr
                << -1 << 4 << 2 << 5 << 5 << endr
                << -2 << 2 << 1 << 1 << -3 << endr;
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_rows); i++){
           for(int j = 0; j < int(sig.n_cols); j++){
            CHECK_EQUAL(real(shifted(i,j)), real(centered(i,j)));
            }
        }
        /*------------------------------------------------------*/





    }

}
