/**********************************************************************
 *  Test: convolution theorem applied on F = W * K = ifft(fft(W)fft(K)),
 *        where W = Gauss(r) * delta(t) and
 *        K(r,t) = Ks(r) * Kt(t) = Cs * Ct
 *
 *  Analytic source: closed-form experssion
 *
 * ********************************************************************/

#include <unittest++/UnitTest++.h>
#include <lgnSimulator.h>

using namespace lgnSimulator;


void runGaussConstConvolutionTest(int nt, double dt, int ns, double ds,
                                int tau, double A, double a,
                                double Cs)
{


    Integrator integrator(nt, dt, ns, ds);
    vec r = integrator.spatialVec();
    vec k = integrator.spatialFreqVec();
    vec t = integrator.timeVec();
    vec w = integrator.temporalFreqVec();

    SpatialGaussian Ws(A, a);
    TemporalDelta Wt(t[tau], dt);

    SpatiallyConstant Ks(Cs, integrator.spatialFreqResolution());
    TemporallyConstant Kt(integrator.timeInterval(),
                          integrator.temporalFreqResolution());


    cube F_e = zeros(r.n_elem, r.n_elem, t.n_elem);
    cx_cube G = zeros<cx_cube>(r.n_elem, r.n_elem, t.n_elem);

    for(int l=0; l < int(t.n_elem); l++){
        for(int i = 0; i < int(r.n_elem); i++){
            for(int j = 0; j < int(r.n_elem); j++){
                F_e(i,j,l) = Cs * (1./integrator.timeInterval()) * A;

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

    TEST(gaussConstConvolutionTest_test_0) {
        runGaussConstConvolutionTest(3, 0.05, 4, 0.05,
                                   2, 1.0, 0.25,
                                   1.0);
    }

    TEST(gaussConstConvolutionTest_test_1) {
        runGaussConstConvolutionTest(2, 0.05, 3, 0.05,
                                   1, -10.263, 0.25,
                                   1.3);
    }


}
