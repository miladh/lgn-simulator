/**********************************************************************
 *  Test: spatial and temporal kernel functions
 *
 *  Analytic source: by hand and Python
 *
 * ********************************************************************/

#include <unittest++/UnitTest++.h>
#include <lgnSimulator.h>

using namespace lgnSimulator;



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

    CHECK_CLOSE(diff_real.max(), 0.0, 1e-7);
    CHECK_CLOSE(diff_imag.max(), 0.0, 1e-7);

}

SUITE(integrator){
    TEST(spatiotemporalDelta_test_0) {
        runSpatiotemporalDeltaTest(6, 0.05, 6, 0.05, 0.0, 0.0,
                                   1.0, 0.25, 1.0, vec2{0.0, 0.0});
    }

    TEST(spatiotemporalDelta_test_1) {
        runSpatiotemporalDeltaTest(7, 0.05, 7, 0.05, 0.3, 1.0,
                                   -1.2, 0.25, 1.0, vec2{1.1, -0.3});
    }


}
