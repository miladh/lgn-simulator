/**********************************************************************
 *  Test: convolution theorem applied on F = W * K = ifft(fft(W)fft(K)),
 *        where W = delta(r) * delta(t) and
 *        K = delta(r) * delta(t)
 *
 *  Analytic source: closed-form experssion
 *
 * ********************************************************************/

#include <unittest++/UnitTest++.h>
#include <lgnSimulator.h>

using namespace lgnSimulator;

void runDeltaConvolutionTest(int nt, double dt, int ns, double ds,
                             int delay_w, vec2 shift_w,
                             int delay_k, vec2 shift_k)
{


    Integrator integrator(nt, dt, ns, ds);
    vec r = integrator.spatialVec();
    vec k = integrator.spatialFreqVec();
    vec t = integrator.timeVec();
    vec w = integrator.temporalFreqVec();

    SpatialDelta Ws(ds, shift_w);
    TemporalDelta Wt(t[delay_w], dt);

    SpatialDelta Ks(ds, shift_k);
    TemporalDelta Kt(t[delay_k], dt);

    cube F_e = zeros(r.n_elem, r.n_elem, t.n_elem);
    cx_cube G = zeros<cx_cube>(r.n_elem, r.n_elem, t.n_elem);


    for(int l=0; l < int(t.n_elem); l++){
        for(int i = 0; i < int(r.n_elem); i++){
            for(int j = 0; j < int(r.n_elem); j++){
                F_e(i,j,l) =  Ws.spatial(vec2{r[i], r[j]} - shift_k)
                           * Wt.temporal(t[l] - t[delay_k]);
                G(i,j,l) = Ws.fourierTransform({k[i], k[j]})
                         * Wt.fourierTransform(w[l])
                         * Ks.fourierTransform({k[i], k[j]})
                         * Kt.fourierTransform(w[l]);

            }
        }
    }

    cube F = integrator.backwardFFT(G);
    cube diff = F - F_e;

    cube diff_real = abs(real(diff));
    cube diff_imag = abs(imag(diff));

    //    cout << diff_real.max() << endl;
    //    cout << diff_imag.max() << endl;

    CHECK_CLOSE(diff_real.max(), 0.0, 1e-10);
    CHECK_CLOSE(diff_imag.max(), 0.0, 1e-10);

}



SUITE(integrator){

    TEST(deltaConvolutionTest_test_0) {
        runDeltaConvolutionTest(2, 0.05, 6, 0.05,
                               0,  vec2{0.0, 0.0},
                               0,  vec2{0.2, 0.0});
    }

    TEST(deltaConvolutionTest_test_1) {
        runDeltaConvolutionTest(3, 0.1, 7, 0.05,
                                0, vec2{0.3, 0.0},
                                0, vec2{0.0, .2});
    }


}
