#include <unittest++/UnitTest++.h>
#include <armadillo>

#include "math/ffthelper.h"


using namespace std;
using namespace arma;
using namespace edog;

SUITE(fftHelper){

    TEST(fftFreq){
        vec result, freqs;

        result = { 0.0, 0.221, 0.441, -0.441, -0.221 };
        freqs = FFTHelper::fftFreq( 5 , 0.907 );
        for(int i = 0; i < int(result.n_elem); i++){
            CHECK_CLOSE(freqs[i], result[i], 1e-3);
        }

        result = { 0.0, 0.187, 0.374, 0.561, -0.749, -0.561, -0.374, -0.187 };
        freqs = FFTHelper::fftFreq( 8 , 0.668 );
        for(int i = 0; i < int(result.n_elem); i++){
            CHECK_CLOSE(freqs[i], result[i], 1e-3);
        }

        result = { 0.0, 0.248, 0.495, -0.495, -0.248 };
        freqs = FFTHelper::fftFreq( 5 , 0.808 );
        for(int i = 0; i < int(result.n_elem); i++){
            CHECK_CLOSE(freqs[i], result[i], 1e-3);
        }

        result = { 0.0, 2.342, 4.684, 7.026, -7.026, -4.684, -2.342 };
        freqs = FFTHelper::fftFreq( 7 , 0.061 );
        for(int i = 0; i < int(result.n_elem); i++){
            CHECK_CLOSE(freqs[i], result[i], 1e-3);
        }

        result = { 0.0, 0.148, 0.295, 0.443, 0.59, -0.59, -0.443, -0.295, -0.148 };
        freqs = FFTHelper::fftFreq( 9 , 0.753 );
        for(int i = 0; i < int(result.n_elem); i++){
            CHECK_CLOSE(freqs[i], result[i], 1e-3);
        }

        result = { 0.0, 0.326, 0.651, -0.651, -0.326 };
        freqs = FFTHelper::fftFreq( 5 , 0.614 );
        for(int i = 0; i < int(result.n_elem); i++){
            CHECK_CLOSE(freqs[i], result[i], 1e-3);
        }

        result = { 0.0, 0.295, 0.589, 0.884, 1.179, -1.179, -0.884, -0.589, -0.295 };
        freqs = FFTHelper::fftFreq( 9 , 0.377 );
        for(int i = 0; i < int(result.n_elem); i++){
            CHECK_CLOSE(freqs[i], result[i], 1e-3);
        }

        result = { 0.0, 0.181, 0.361, 0.542, 0.723, -0.723, -0.542, -0.361, -0.181 };
        freqs = FFTHelper::fftFreq( 9 , 0.615 );
        for(int i = 0; i < int(result.n_elem); i++){
            CHECK_CLOSE(freqs[i], result[i], 1e-3);
        }

        result = { 0.0, 1.106, 2.212, 3.319, -4.425, -3.319, -2.212, -1.106 };
        freqs = FFTHelper::fftFreq( 8 , 0.113 );
        for(int i = 0; i < int(result.n_elem); i++){
            CHECK_CLOSE(freqs[i], result[i], 1e-3);
        }

        result = { 0.0, 0.329, 0.659, -0.659, -0.329 };
        freqs = FFTHelper::fftFreq( 5 , 0.607 );
        for(int i = 0; i < int(result.n_elem); i++){
            CHECK_CLOSE(freqs[i], result[i], 1e-3);
        }

        result = { 0.0, 0.443, 0.885, 1.328, 1.771, -1.771, -1.328, -0.885, -0.443 };
        freqs = FFTHelper::fftFreq( 9 , 0.251 );
        for(int i = 0; i < int(result.n_elem); i++){
            CHECK_CLOSE(freqs[i], result[i], 1e-3);
        }

        result = { 0.0, 0.567, 1.134, -1.701, -1.134, -0.567 };
        freqs = FFTHelper::fftFreq( 6 , 0.294 );
        for(int i = 0; i < int(result.n_elem); i++){
            CHECK_CLOSE(freqs[i], result[i], 1e-3);
        }

        result = { 0.0, 0.342, 0.685, 1.027, -1.37, -1.027, -0.685, -0.342 };
        freqs = FFTHelper::fftFreq( 8 , 0.365 );
        for(int i = 0; i < int(result.n_elem); i++){
            CHECK_CLOSE(freqs[i], result[i], 1e-3);
        }

        result = { 0.0, 0.14, 0.279, 0.419, 0.558, -0.558, -0.419, -0.279, -0.14 };
        freqs = FFTHelper::fftFreq( 9 , 0.796 );
        for(int i = 0; i < int(result.n_elem); i++){
            CHECK_CLOSE(freqs[i], result[i], 1e-3);
        }

        result = { 0.0, 0.39, 0.78, 1.17, 1.559, -1.559, -1.17, -0.78, -0.39 };
        freqs = FFTHelper::fftFreq( 9 , 0.285 );
        for(int i = 0; i < int(result.n_elem); i++){
            CHECK_CLOSE(freqs[i], result[i], 1e-3);
        }

        result = { 0.0, 0.621, 1.242, 1.863, -1.863, -1.242, -0.621 };
        freqs = FFTHelper::fftFreq( 7 , 0.23 );
        for(int i = 0; i < int(result.n_elem); i++){
            CHECK_CLOSE(freqs[i], result[i], 1e-3);
        }

        result = { 0.0, 0.329, 0.657, -0.986, -0.657, -0.329 };
        freqs = FFTHelper::fftFreq( 6 , 0.507 );
        for(int i = 0; i < int(result.n_elem); i++){
            CHECK_CLOSE(freqs[i], result[i], 1e-3);
        }

        result = { 0.0, 0.118, 0.236, 0.354, 0.472, -0.472, -0.354, -0.236, -0.118 };
        freqs = FFTHelper::fftFreq( 9 , 0.941 );
        for(int i = 0; i < int(result.n_elem); i++){
            CHECK_CLOSE(freqs[i], result[i], 1e-3);
        }

        result = { 0.0, 2.165, 4.329, -6.494, -4.329, -2.165 };
        freqs = FFTHelper::fftFreq( 6 , 0.077 );
        for(int i = 0; i < int(result.n_elem); i++){
            CHECK_CLOSE(freqs[i], result[i], 1e-3);
        }

        result = { 0.0, 0.176, 0.352, 0.528, -0.528, -0.352, -0.176 };
        freqs = FFTHelper::fftFreq( 7 , 0.812 );
        for(int i = 0; i < int(result.n_elem); i++){
            CHECK_CLOSE(freqs[i], result[i], 1e-3);
        }

        result = { 0.0, 10.204, 20.408, 30.612, -30.612, -20.408, -10.204 };
        freqs = FFTHelper::fftFreq( 7 , 0.014 );
        for(int i = 0; i < int(result.n_elem); i++){
            CHECK_CLOSE(freqs[i], result[i], 1e-3);
        }

        result = { 0.0, 0.505, 1.01, -1.01, -0.505 };
        freqs = FFTHelper::fftFreq( 5 , 0.396 );
        for(int i = 0; i < int(result.n_elem); i++){
            CHECK_CLOSE(freqs[i], result[i], 1e-3);
        }

        result = { 0.0, 0.358, 0.717, 1.075, 1.434, -1.434, -1.075, -0.717, -0.358 };
        freqs = FFTHelper::fftFreq( 9 , 0.31 );
        for(int i = 0; i < int(result.n_elem); i++){
            CHECK_CLOSE(freqs[i], result[i], 1e-3);
        }

        result = { 0.0, 0.207, 0.415, 0.622, -0.622, -0.415, -0.207 };
        freqs = FFTHelper::fftFreq( 7 , 0.689 );
        for(int i = 0; i < int(result.n_elem); i++){
            CHECK_CLOSE(freqs[i], result[i], 1e-3);
        }

        result = { 0.0, 14.286, 28.571, 42.857, -42.857, -28.571, -14.286 };
        freqs = FFTHelper::fftFreq( 7 , 0.01 );
        for(int i = 0; i < int(result.n_elem); i++){
            CHECK_CLOSE(freqs[i], result[i], 1e-3);
        }

        result = { 0.0, 0.225, 0.45, 0.675, 0.9, -0.9, -0.675, -0.45, -0.225 };
        freqs = FFTHelper::fftFreq( 9 , 0.494 );
        for(int i = 0; i < int(result.n_elem); i++){
            CHECK_CLOSE(freqs[i], result[i], 1e-3);
        }

        result = { 0.0, 0.229, 0.457, 0.686, 0.914, -0.914, -0.686, -0.457, -0.229 };
        freqs = FFTHelper::fftFreq( 9 , 0.486 );
        for(int i = 0; i < int(result.n_elem); i++){
            CHECK_CLOSE(freqs[i], result[i], 1e-3);
        }

        result = { 0.0, 0.118, 0.235, 0.353, 0.47, -0.47, -0.353, -0.235, -0.118 };
        freqs = FFTHelper::fftFreq( 9 , 0.945 );
        for(int i = 0; i < int(result.n_elem); i++){
            CHECK_CLOSE(freqs[i], result[i], 1e-3);
        }

        result = { 0.0, 0.121, 0.241, 0.362, 0.483, -0.483, -0.362, -0.241, -0.121 };
        freqs = FFTHelper::fftFreq( 9 , 0.921 );
        for(int i = 0; i < int(result.n_elem); i++){
            CHECK_CLOSE(freqs[i], result[i], 1e-3);
        }

        result = { 0.0, 0.212, 0.424, -0.635, -0.424, -0.212 };
        freqs = FFTHelper::fftFreq( 6 , 0.787 );
        for(int i = 0; i < int(result.n_elem); i++){
            CHECK_CLOSE(freqs[i], result[i], 1e-3);
        }




    }

}
