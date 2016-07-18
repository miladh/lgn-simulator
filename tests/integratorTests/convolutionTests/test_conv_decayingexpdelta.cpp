/**********************************************************************
 *  Test: convolution theorem applied on F = W * K = ifft(fft(W)fft(K)),
 *        where W = Gauss(r) * delta(t) and
 *        K = delta(r) * DecayingExponential(t)
 *
 *  Analytic source: closed-form experssion
 *
 * ********************************************************************/

#include <unittest++/UnitTest++.h>
#include <lgnSimulator.h>

using namespace lgnSimulator;


void runDogdecayingExpConvolutionTest(int nt, double dt, int ns, double ds,
                                      int tau, double a,
                                      int delay, double timeConstant,
                                      vec2 shift)
{


    Integrator integrator(nt, dt, ns, ds);
    vec r = integrator.spatialVec();
    vec k = integrator.spatialFreqVec();
    vec t = integrator.timeVec();
    vec w = integrator.temporalFreqVec();

    SpatialGaussian Ws( a);
    TemporalDelta Wt(t[tau], dt);

    SpatialDelta Ks( ds, shift);
    DecayingExponential Kt(timeConstant, t[delay]);

    cube F_e = zeros(r.n_elem, r.n_elem, t.n_elem);
    cx_cube G = zeros<cx_cube>(r.n_elem, r.n_elem, t.n_elem);

    for(int l=0; l < int(t.n_elem); l++){
        for(int i = 0; i < int(r.n_elem); i++){
            for(int j = 0; j < int(r.n_elem); j++){
                F_e(i,j,l) =  Ws.spatial(vec2{r[i], r[j]} - shift)
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

//    cout << diff_real.max() << endl;
//    cout << diff_imag.max() << endl;

    CHECK_CLOSE(diff_real.max(), 0.0, 1e-7);
    CHECK_CLOSE(diff_imag.max(), 0.0, 1e-7);

}


SUITE(integrator){

//    TEST(runDogdecayingExpConvolutionTest_test_0) {
//        runDogdecayingExpConvolutionTest(10, 1.0, 6, 0.05,
//                                         0, 1.0, 0.25,
//                                         0, 0.04, 1.0,
//                                         vec2{0.0, 0.0});
//    }

}
