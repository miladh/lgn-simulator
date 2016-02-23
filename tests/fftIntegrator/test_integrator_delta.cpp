/**********************************************************************
 *  Test: spatial and temporal kernel functions
 *
 *  Analytic source: by hand and Python
 *
 * ********************************************************************/

#include <unittest++/UnitTest++.h>
#include <lgnSimulator.h>

using namespace lgnSimulator;


void runTemporalDeltaTest(int ns, double ds,
                           double A, double a, double w, vec2 shift){


    Gaussian G(A, a);
    SpatialDelta delta(w, shift);

    Integrator integrator(2, 0.1, ns, ds);
    vec r = integrator.spatialVec();
    vec k = integrator.spatialFreqVec();
    mat spatial = zeros(r.n_elem, r.n_elem);
    cx_mat fourierTransform = zeros<cx_mat>(r.n_elem, r.n_elem);

    for(int i = 0; i < int(r.n_elem); i++){
        for(int j = 0; j < int(r.n_elem); j++){
            spatial(i,j) = w * G.spatial(vec2{r[i], r[j]}-shift);
            fourierTransform(i,j) = G.fourierTransform({k[i], k[j]})
                    * delta.fourierTransform({k[i], k[j]});
        }
    }

    cx_mat spatial_fft = integrator.backwardFFT(fourierTransform);
    cx_mat diff = spatial_fft - spatial;

    mat diff_real = abs(real(diff));
    mat diff_imag = abs(imag(diff));

    CHECK_CLOSE(diff_real.max(), 0.0, 1e-9);
    CHECK_CLOSE(diff_imag.max(), 0.0, 1e-9);

}



SUITE(integrator){
    TEST(temporalDelta_test_0) {
        runTemporalDeltaTest(7, 0.05,  1.0, 0.25, 1.0, vec2{0.0, 0.0});
    }

    TEST(temporalDelta_test_1) {
        runTemporalDeltaTest(7, 0.05,  1.0, 0.35, 2.4, vec2{0.2, 1.4});
    }

    TEST(temporalDelta_test_2) {
        runTemporalDeltaTest(9, 0.05,  -1.2, -0.65, -1.1, vec2{-1.2, 8.4});
    }

    TEST(temporalDelta_test_3) {
        runTemporalDeltaTest(9, 0.01,  4.2, 0.07, 4.2, vec2{0.6, -1.4});
    }

}
