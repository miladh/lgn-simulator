
#include <unittest++/UnitTest++.h>
#include <armadillo>

#include "math/ffthelper.h"


using namespace std;
using namespace arma;

SUITE(fftHelper){


    TEST(fftshift_1D){
        cx_vec sig, centered, shifted;


        sig = { 0.0, 0.184, 0.368, 0.552, 0.736, -0.736, -0.552, -0.368, -0.184 };
        centered = { -0.736, -0.552, -0.368, -0.184, 0.0, 0.184, 0.368, 0.552, 0.736 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_CLOSE(real(shifted[i]), real(centered[i]), 1e-3);
        }

        sig = { 0.0, 0.272, 0.545, 0.817, 1.089, -1.089, -0.817, -0.545, -0.272 };
        centered = { -1.089, -0.817, -0.545, -0.272, 0.0, 0.272, 0.545, 0.817, 1.089 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_CLOSE(real(shifted[i]), real(centered[i]), 1e-3);
        }

        sig = { 0.0, 0.324, 0.648, 0.972, 1.296, -1.296, -0.972, -0.648, -0.324 };
        centered = { -1.296, -0.972, -0.648, -0.324, 0.0, 0.324, 0.648, 0.972, 1.296 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_CLOSE(real(shifted[i]), real(centered[i]), 1e-3);
        }

        sig = { 0.0, 0.427, 0.855, 1.282, 1.709, -1.709, -1.282, -0.855, -0.427 };
        centered = { -1.709, -1.282, -0.855, -0.427, 0.0, 0.427, 0.855, 1.282, 1.709 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_CLOSE(real(shifted[i]), real(centered[i]), 1e-3);
        }

        sig = { 0.0, 0.197, 0.395, 0.592, -0.592, -0.395, -0.197 };
        centered = { -0.592, -0.395, -0.197, 0.0, 0.197, 0.395, 0.592 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_CLOSE(real(shifted[i]), real(centered[i]), 1e-3);
        }

        sig = { 0.0, 0.361, 0.723, 1.084, -1.445, -1.084, -0.723, -0.361 };
        centered = { -1.445, -1.084, -0.723, -0.361, 0.0, 0.361, 0.723, 1.084 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_CLOSE(real(shifted[i]), real(centered[i]), 1e-3);
        }

        sig = { 0.0, 0.174, 0.348, 0.522, 0.697, -0.697, -0.522, -0.348, -0.174 };
        centered = { -0.697, -0.522, -0.348, -0.174, 0.0, 0.174, 0.348, 0.522, 0.697 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_CLOSE(real(shifted[i]), real(centered[i]), 1e-3);
        }

        sig = { 0.0, 0.604, 1.208, -1.208, -0.604 };
        centered = { -1.208, -0.604, 0.0, 0.604, 1.208 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_CLOSE(real(shifted[i]), real(centered[i]), 1e-3);
        }

        sig = { 0.0, 0.195, 0.39, 0.585, -0.585, -0.39, -0.195 };
        centered = { -0.585, -0.39, -0.195, 0.0, 0.195, 0.39, 0.585 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_CLOSE(real(shifted[i]), real(centered[i]), 1e-3);
        }

        sig = { 0.0, 0.289, 0.577, -0.577, -0.289 };
        centered = { -0.577, -0.289, 0.0, 0.289, 0.577 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_CLOSE(real(shifted[i]), real(centered[i]), 1e-3);
        }

        sig = { 0.0, 1.124, 2.247, -2.247, -1.124 };
        centered = { -2.247, -1.124, 0.0, 1.124, 2.247 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_CLOSE(real(shifted[i]), real(centered[i]), 1e-3);
        }

        sig = { 0.0, 0.234, 0.468, -0.702, -0.468, -0.234 };
        centered = { -0.702, -0.468, -0.234, 0.0, 0.234, 0.468 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_CLOSE(real(shifted[i]), real(centered[i]), 1e-3);
        }

        sig = { 0.0, 0.152, 0.304, 0.456, -0.456, -0.304, -0.152 };
        centered = { -0.456, -0.304, -0.152, 0.0, 0.152, 0.304, 0.456 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_CLOSE(real(shifted[i]), real(centered[i]), 1e-3);
        }

        sig = { 0.0, 0.699, 1.399, -1.399, -0.699 };
        centered = { -1.399, -0.699, 0.0, 0.699, 1.399 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_CLOSE(real(shifted[i]), real(centered[i]), 1e-3);
        }

        sig = { 0.0, 0.465, 0.929, 1.394, -1.859, -1.394, -0.929, -0.465 };
        centered = { -1.859, -1.394, -0.929, -0.465, 0.0, 0.465, 0.929, 1.394 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_CLOSE(real(shifted[i]), real(centered[i]), 1e-3);
        }

        sig = { 0.0, 0.362, 0.723, -0.723, -0.362 };
        centered = { -0.723, -0.362, 0.0, 0.362, 0.723 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_CLOSE(real(shifted[i]), real(centered[i]), 1e-3);
        }

        sig = { 0.0, 0.244, 0.488, -0.732, -0.488, -0.244 };
        centered = { -0.732, -0.488, -0.244, 0.0, 0.244, 0.488 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_CLOSE(real(shifted[i]), real(centered[i]), 1e-3);
        }

        sig = { 0.0, 0.116, 0.232, 0.349, 0.465, -0.465, -0.349, -0.232, -0.116 };
        centered = { -0.465, -0.349, -0.232, -0.116, 0.0, 0.116, 0.232, 0.349, 0.465 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_CLOSE(real(shifted[i]), real(centered[i]), 1e-3);
        }

        sig = { 0.0, 5.556, 11.111, 16.667, 22.222, -22.222, -16.667, -11.111, -5.556 };
        centered = { -22.222, -16.667, -11.111, -5.556, 0.0, 5.556, 11.111, 16.667, 22.222 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_CLOSE(real(shifted[i]), real(centered[i]), 1e-3);
        }

        sig = { 0.0, 0.192, 0.383, 0.575, -0.767, -0.575, -0.383, -0.192 };
        centered = { -0.767, -0.575, -0.383, -0.192, 0.0, 0.192, 0.383, 0.575 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_CLOSE(real(shifted[i]), real(centered[i]), 1e-3);
        }

        sig = { 0.0, 0.458, 0.916, -1.374, -0.916, -0.458 };
        centered = { -1.374, -0.916, -0.458, 0.0, 0.458, 0.916 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_CLOSE(real(shifted[i]), real(centered[i]), 1e-3);
        }

        sig = { 0.0, 0.676, 1.351, -1.351, -0.676 };
        centered = { -1.351, -0.676, 0.0, 0.676, 1.351 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_CLOSE(real(shifted[i]), real(centered[i]), 1e-3);
        }

        sig = { 0.0, 0.195, 0.39, -0.585, -0.39, -0.195 };
        centered = { -0.585, -0.39, -0.195, 0.0, 0.195, 0.39 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_CLOSE(real(shifted[i]), real(centered[i]), 1e-3);
        }

        sig = { 0.0, 0.182, 0.365, -0.547, -0.365, -0.182 };
        centered = { -0.547, -0.365, -0.182, 0.0, 0.182, 0.365 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_CLOSE(real(shifted[i]), real(centered[i]), 1e-3);
        }

        sig = { 0.0, 0.252, 0.504, -0.504, -0.252 };
        centered = { -0.504, -0.252, 0.0, 0.252, 0.504 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_CLOSE(real(shifted[i]), real(centered[i]), 1e-3);
        }

        sig = { 0.0, 0.128, 0.256, 0.384, 0.513, -0.513, -0.384, -0.256, -0.128 };
        centered = { -0.513, -0.384, -0.256, -0.128, 0.0, 0.128, 0.256, 0.384, 0.513 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_CLOSE(real(shifted[i]), real(centered[i]), 1e-3);
        }

        sig = { 0.0, 0.163, 0.326, 0.489, -0.489, -0.326, -0.163 };
        centered = { -0.489, -0.326, -0.163, 0.0, 0.163, 0.326, 0.489 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_CLOSE(real(shifted[i]), real(centered[i]), 1e-3);
        }

        sig = { 0.0, 6.944, 13.889, -20.833, -13.889, -6.944 };
        centered = { -20.833, -13.889, -6.944, 0.0, 6.944, 13.889 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_CLOSE(real(shifted[i]), real(centered[i]), 1e-3);
        }

        sig = { 0.0, 0.291, 0.583, -0.874, -0.583, -0.291 };
        centered = { -0.874, -0.583, -0.291, 0.0, 0.291, 0.583 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_CLOSE(real(shifted[i]), real(centered[i]), 1e-3);
        }

        sig = { 0.0, 0.329, 0.658, 0.987, -0.987, -0.658, -0.329 };
        centered = { -0.987, -0.658, -0.329, 0.0, 0.329, 0.658, 0.987 };
        shifted = FFTHelper::fftShift(sig);
        for(int i = 0; i < int(sig.n_elem); i++){
            CHECK_CLOSE(real(shifted[i]), real(centered[i]), 1e-3);
        }



    }

}
