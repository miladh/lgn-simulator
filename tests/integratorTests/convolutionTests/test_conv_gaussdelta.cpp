/**********************************************************************
 *  Test: convolution theorem applied on F = W * K = ifft(fft(W)fft(K)),
 *        where W = Gauss(r) * delta(t) and
 *        K = delta(r) * delta(t)
 *
 *  Analytic source: closed-form experssion
 *
 * ********************************************************************/


#include <lgnSimulator.h>
#include <catch.hpp>


using namespace lgnSimulator;


void runGaussDeltaConvolutionTest(int nt, double dt, int ns, double ds,
                                int tau, double a,
                                int delay,  vec2 shift)
{


    Integrator integrator(nt, dt, ns, ds);
    vec r = integrator.spatialVec();
    vec k = integrator.spatialFreqVec();
    vec t = integrator.timeVec();
    vec w = integrator.temporalFreqVec();

    SpatialGaussian Ws( a);
    TemporalDelta Wt(t[tau], dt);

    SpatialDelta Ks(ds, shift);
    TemporalDelta Kt(t[delay], dt);

    cube F_e = zeros(r.n_elem, r.n_elem, t.n_elem);
    cx_cube G = zeros<cx_cube>(r.n_elem, r.n_elem, t.n_elem);


    for(int l=0; l < int(t.n_elem); l++){
        for(int i = 0; i < int(r.n_elem); i++){
            for(int j = 0; j < int(r.n_elem); j++){
                F_e(i,j,l) = Ws.spatial(vec2{r[i], r[j]} - shift)
                           * Wt.temporal(t[l] - t[delay]);
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

    REQUIRE(diff_real.max()==  Approx( 0.0).epsilon(1e-7));
    REQUIRE(diff_imag.max()==  Approx( 0.0).epsilon(1e-7));

}


TEST_CASE("gaussDeltaConvolutionTest_test_0") {
    runGaussDeltaConvolutionTest(2, 0.05, 6, 0.05,
                                 0,  0.25,
                                 0,  vec2{0.0, 0.0});
}

TEST_CASE("gaussDeltaConvolutionTest_test_1") {
    runGaussDeltaConvolutionTest(3, 0.1, 7, 0.05,
                                 1,  0.25,
                                 5,  vec2{1.1, -0.3});
}

