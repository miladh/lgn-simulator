
#include <unittest++/UnitTest++.h>
#include <armadillo>

#include "math/ffthelper.h"


using namespace std;
using namespace arma;
using namespace edog;

SUITE(fftHelper){


    TEST(fftshift_1D){
        cx_vec sig, centered, shifted;


        sig = { 0.0, 1.46, 2.92, -2.92, -1.46 };
        centered = { -2.92, -1.46, 0.0, 1.46, 2.92 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_EQUAL(real(shifted[i]), real(centered[i]));
        }

        sig = { 0.0, 0.121, 0.242, 0.362, 0.483, -0.483, -0.362, -0.242, -0.121 };
        centered = { -0.483, -0.362, -0.242, -0.121, 0.0, 0.121, 0.242, 0.362, 0.483 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_EQUAL(real(shifted[i]), real(centered[i]));
        }

        sig = { 0.0, 0.483, 0.965, 1.448, -1.448, -0.965, -0.483 };
        centered = { -1.448, -0.965, -0.483, 0.0, 0.483, 0.965, 1.448 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_EQUAL(real(shifted[i]), real(centered[i]));
        }

        sig = { 0.0, 1.068, 2.137, 3.205, 4.274, -4.274, -3.205, -2.137, -1.068 };
        centered = { -4.274, -3.205, -2.137, -1.068, 0.0, 1.068, 2.137, 3.205, 4.274 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_EQUAL(real(shifted[i]), real(centered[i]));
        }

        sig = { 0.0, 0.151, 0.303, 0.454, -0.605, -0.454, -0.303, -0.151 };
        centered = { -0.605, -0.454, -0.303, -0.151, 0.0, 0.151, 0.303, 0.454 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_EQUAL(real(shifted[i]), real(centered[i]));
        }

        sig = { 0.0, 1.949, 3.899, 5.848, 7.797, -7.797, -5.848, -3.899, -1.949 };
        centered = { -7.797, -5.848, -3.899, -1.949, 0.0, 1.949, 3.899, 5.848, 7.797 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_EQUAL(real(shifted[i]), real(centered[i]));
        }

        sig = { 0.0, 0.239, 0.477, -0.477, -0.239 };
        centered = { -0.477, -0.239, 0.0, 0.239, 0.477 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_EQUAL(real(shifted[i]), real(centered[i]));
        }

        sig = { 0.0, 0.166, 0.332, 0.498, 0.663, -0.663, -0.498, -0.332, -0.166 };
        centered = { -0.663, -0.498, -0.332, -0.166, 0.0, 0.166, 0.332, 0.498, 0.663 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_EQUAL(real(shifted[i]), real(centered[i]));
        }

        sig = { 0.0, 0.253, 0.506, -0.506, -0.253 };
        centered = { -0.506, -0.253, 0.0, 0.253, 0.506 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_EQUAL(real(shifted[i]), real(centered[i]));
        }

        sig = { 0.0, 0.164, 0.327, 0.491, -0.654, -0.491, -0.327, -0.164 };
        centered = { -0.654, -0.491, -0.327, -0.164, 0.0, 0.164, 0.327, 0.491 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_EQUAL(real(shifted[i]), real(centered[i]));
        }

        sig = { 0.0, 0.243, 0.486, 0.73, -0.973, -0.73, -0.486, -0.243 };
        centered = { -0.973, -0.73, -0.486, -0.243, 0.0, 0.243, 0.486, 0.73 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_EQUAL(real(shifted[i]), real(centered[i]));
        }

        sig = { 0.0, 0.144, 0.289, 0.433, -0.433, -0.289, -0.144 };
        centered = { -0.433, -0.289, -0.144, 0.0, 0.144, 0.289, 0.433 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_EQUAL(real(shifted[i]), real(centered[i]));
        }

        sig = { 0.0, 0.127, 0.254, 0.381, 0.508, -0.508, -0.381, -0.254, -0.127 };
        centered = { -0.508, -0.381, -0.254, -0.127, 0.0, 0.127, 0.254, 0.381, 0.508 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_EQUAL(real(shifted[i]), real(centered[i]));
        }

        sig = { 0.0, 0.135, 0.271, 0.406, -0.541, -0.406, -0.271, -0.135 };
        centered = { -0.541, -0.406, -0.271, -0.135, 0.0, 0.135, 0.271, 0.406 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_EQUAL(real(shifted[i]), real(centered[i]));
        }

        sig = { 0.0, 0.593, 1.186, 1.778, -1.778, -1.186, -0.593 };
        centered = { -1.778, -1.186, -0.593, 0.0, 0.593, 1.186, 1.778 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_EQUAL(real(shifted[i]), real(centered[i]));
        }

        sig = { 0.0, 0.131, 0.261, 0.392, 0.522, -0.522, -0.392, -0.261, -0.131 };
        centered = { -0.522, -0.392, -0.261, -0.131, 0.0, 0.131, 0.261, 0.392, 0.522 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_EQUAL(real(shifted[i]), real(centered[i]));
        }

        sig = { 0.0, 0.445, 0.89, 1.335, -1.335, -0.89, -0.445 };
        centered = { -1.335, -0.89, -0.445, 0.0, 0.445, 0.89, 1.335 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_EQUAL(real(shifted[i]), real(centered[i]));
        }

        sig = { 0.0, 0.241, 0.483, 0.724, -0.965, -0.724, -0.483, -0.241 };
        centered = { -0.965, -0.724, -0.483, -0.241, 0.0, 0.241, 0.483, 0.724 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_EQUAL(real(shifted[i]), real(centered[i]));
        }

        sig = { 0.0, 0.146, 0.292, 0.437, -0.437, -0.292, -0.146 };
        centered = { -0.437, -0.292, -0.146, 0.0, 0.146, 0.292, 0.437 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_EQUAL(real(shifted[i]), real(centered[i]));
        }

        sig = { 0.0, 0.226, 0.451, 0.677, -0.677, -0.451, -0.226 };
        centered = { -0.677, -0.451, -0.226, 0.0, 0.226, 0.451, 0.677 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_EQUAL(real(shifted[i]), real(centered[i]));
        }

        sig = { 0.0, 0.386, 0.772, 1.157, -1.543, -1.157, -0.772, -0.386 };
        centered = { -1.543, -1.157, -0.772, -0.386, 0.0, 0.386, 0.772, 1.157 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_EQUAL(real(shifted[i]), real(centered[i]));
        }

        sig = { 0.0, 0.526, 1.053, -1.053, -0.526 };
        centered = { -1.053, -0.526, 0.0, 0.526, 1.053 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_EQUAL(real(shifted[i]), real(centered[i]));
        }

        sig = { 0.0, 0.868, 1.736, -2.604, -1.736, -0.868 };
        centered = { -2.604, -1.736, -0.868, 0.0, 0.868, 1.736 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_EQUAL(real(shifted[i]), real(centered[i]));
        }

        sig = { 0.0, 0.41, 0.82, 1.23, 1.64, -1.64, -1.23, -0.82, -0.41 };
        centered = { -1.64, -1.23, -0.82, -0.41, 0.0, 0.41, 0.82, 1.23, 1.64 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_EQUAL(real(shifted[i]), real(centered[i]));
        }

        sig = { 0.0, 1.244, 2.488, -3.731, -2.488, -1.244 };
        centered = { -3.731, -2.488, -1.244, 0.0, 1.244, 2.488 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_EQUAL(real(shifted[i]), real(centered[i]));
        }

        sig = { 0.0, 0.533, 1.066, 1.599, -1.599, -1.066, -0.533 };
        centered = { -1.599, -1.066, -0.533, 0.0, 0.533, 1.066, 1.599 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_EQUAL(real(shifted[i]), real(centered[i]));
        }

        sig = { 0.0, 0.148, 0.297, 0.445, -0.594, -0.445, -0.297, -0.148 };
        centered = { -0.594, -0.445, -0.297, -0.148, 0.0, 0.148, 0.297, 0.445 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_EQUAL(real(shifted[i]), real(centered[i]));
        }

        sig = { 0.0, 1.157, 2.315, 3.472, 4.63, -4.63, -3.472, -2.315, -1.157 };
        centered = { -4.63, -3.472, -2.315, -1.157, 0.0, 1.157, 2.315, 3.472, 4.63 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_EQUAL(real(shifted[i]), real(centered[i]));
        }

        sig = { 0.0, 0.269, 0.538, 0.807, -0.807, -0.538, -0.269 };
        centered = { -0.807, -0.538, -0.269, 0.0, 0.269, 0.538, 0.807 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_EQUAL(real(shifted[i]), real(centered[i]));
        }

        sig = { 0.0, 0.185, 0.37, 0.555, -0.74, -0.555, -0.37, -0.185 };
        centered = { -0.74, -0.555, -0.37, -0.185, 0.0, 0.185, 0.37, 0.555 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_EQUAL(real(shifted[i]), real(centered[i]));
        }


    }

}
