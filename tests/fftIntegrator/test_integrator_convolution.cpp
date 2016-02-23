/**********************************************************************
 *  Test: spatial and temporal kernel functions
 *
 *  Analytic source: by hand and Python
 *
 * ********************************************************************/

#include <unittest++/UnitTest++.h>
#include <lgnSimulator.h>

using namespace lgnSimulator;





void runSpatialDeltaTest(int ns, double ds,
                         double A, double a, double w, vec2 shift){


    Gaussian G(A, a);
    SpatialDelta delta(w, ds, shift);

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

void runSpatiotemporalDeltaTest(int nt, double dt, int ns, double ds,
                                int tau, int delay,
                                double A, double a, double weight, vec2 shift
                                ){


    Integrator integrator(nt, dt, ns, ds);
    vec r = integrator.spatialVec();
    vec k = integrator.spatialFreqVec();
    vec t = integrator.timeVec();
    vec w = integrator.temporalFreqVec();

    Gaussian Ws(A, a);
    TemporalDelta Wt(t[tau], dt);

    SpatialDelta Ks(weight, ds, shift);
    TemporalDelta Kt(t[delay], dt);


    cube F_e = zeros(r.n_elem, r.n_elem, t.n_elem);
    cx_cube G = zeros<cx_cube>(r.n_elem, r.n_elem, t.n_elem);

    for(int l=0; l < int(k.n_elem); l++){
        for(int i = 0; i < int(r.n_elem); i++){
            for(int j = 0; j < int(r.n_elem); j++){
                F_e(i,j,l) = weight
                           * Ws.spatial(vec2{r[i], r[j]} - shift)
                           * Wt.temporal(t[l] - t[delay]);

                G(i,j,l) = Ws.fourierTransform({k[i], k[j]})
                         * Wt.fourierTransform(w[l])
                         * Ks.fourierTransform({k[i], k[j]})
                         * Kt.fourierTransform(w[l]);
            }
        }
    }


    cx_cube F = integrator.backwardFFT(G);
    cx_cube diff = F - F_e;

    cube diff_real = abs(real(diff));
    cube diff_imag = abs(imag(diff));

    cout << diff_real.max() << endl;
    cout << diff_imag.max() << endl << endl;


    CHECK_CLOSE(diff_real.max(), 0.0, 1e-7);
    CHECK_CLOSE(diff_imag.max(), 0.0, 1e-7);

}

SUITE(integrator){
    TEST(spatiotemporalDelta_test_0) {
        runSpatiotemporalDeltaTest(7, 0.05, 7, 0.05, 0.0, 0.0,
                                   1.0, 0.25, 1.0, vec2{0.0, 0.0});
    }

    TEST(spatiotemporalDelta_test_1) {
        runSpatiotemporalDeltaTest(7, 0.05, 7, 0.05, 0.3, 1.0,
                                   1.0, 0.25, 1.0, vec2{0.0, 0.0});
    }

    TEST(spatialDelta_test_0) {
        runSpatialDeltaTest(7, 0.05,  1.0, 0.25, 1.0, vec2{0.0, 0.0});
    }

    TEST(spatialDelta_test_1) {
        runSpatialDeltaTest(7, 0.05,  1.0, 0.35, 2.4, vec2{0.2, 1.4});
    }

    TEST(spatialDelta_test_2) {
        runSpatialDeltaTest(9, 0.05,  -1.2, -0.65, -1.1, vec2{-1.2, 9.8});
    }

    TEST(spatialDelta_test_3) {
        runSpatialDeltaTest(9, 0.01,  4.2, 0.07, 4.2, vec2{0.6, -1.4});
    }

}
