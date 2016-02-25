/**********************************************************************
 *  Test: convolution theorem applied on F = W * K = ifft(fft(W)fft(K)),
 *        where W = Gauss(r) * delta(t) and
 *        K = delta(r) * delta(t)
 *
 *  Analytic source: closed-form experssion
 *
 * ********************************************************************/

#include <unittest++/UnitTest++.h>
#include <lgnSimulator.h>

using namespace lgnSimulator;

void runDogDeltaConvolutionTest(int nt, double dt, int ns, double ds,
                                int tau, double A, double a,
                                int delay, double weight, vec2 shift)
{


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


    for(int l=0; l < int(t.n_elem); l++){
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


void runDogConstConvolutionTest(int nt, double dt, int ns, double ds,
                                int tau, double A, double a,
                                double Cs, double Ct)
{


    Integrator integrator(nt, dt, ns, ds);
    vec r = integrator.spatialVec();
    vec k = integrator.spatialFreqVec();
    vec t = integrator.timeVec();
    vec w = integrator.temporalFreqVec();

    Gaussian Ws(A, a);
    TemporalDelta Wt(t[tau], dt);

    SpatiallyConstant Ks(Cs, integrator.spatialFreqResolution());
    TemporallyConstant Kt(Ct, integrator.temporalFreqResolution());


    cube F_e = zeros(r.n_elem, r.n_elem, t.n_elem);
    cx_cube G = zeros<cx_cube>(r.n_elem, r.n_elem, t.n_elem);

    for(int l=0; l < int(t.n_elem); l++){
        for(int i = 0; i < int(r.n_elem); i++){
            for(int j = 0; j < int(r.n_elem); j++){
                F_e(i,j,l) = Cs * Ct * A;

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

    TEST(dogConstConvolutionTest_test_0) {
        runDogConstConvolutionTest(3, 0.05, 4, 0.05,
                                   2, 1.0, 0.25,
                                   1.0, 2.0);
    }

    TEST(dogConstConvolutionTest_test_1) {
        runDogConstConvolutionTest(2, 0.05, 3, 0.05,
                                   1, -10.263, 0.25,
                                   1.3, -2000.467);
    }

    TEST(dogDeltaConvolutionTest_test_0) {
        runDogDeltaConvolutionTest(2, 0.05, 6, 0.05,
                                   0, 1.0, 0.25,
                                   0, 1.0, vec2{0.0, 0.0});
    }

    TEST(dogDeltaConvolutionTest_test_1) {
        runDogDeltaConvolutionTest(3, 0.1, 7, 0.05,
                                   1, -1.2, 0.25,
                                   5, 10.89, vec2{1.1, -0.3});
    }



}
