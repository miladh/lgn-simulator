/**********************************************************************
 *  Test: convolution theorem applied on F = W * K = ifft(fft(W)fft(K)),
 *        where W = Gauss(r) * delta(t) and
 *        K(r,t) = Ks(r) * Kt(t) = Cs * Ct
 *
 *  Analytic source: closed-form experssion
 *
 * ********************************************************************/


#include <lgnSimulator.h>
#include <catch.hpp>


using namespace lgnSimulator;



void runGaussConstConvolutionTest(int nt, double dt, int ns, double ds,
                                  int tau, double a)
{


    Integrator integrator(nt, dt, ns, ds);
    vec r = integrator.spatialVec();
    vec k = integrator.spatialFreqVec();
    vec t = integrator.timeVec();
    vec w = integrator.temporalFreqVec();

    SpatialGaussian Ws(a);
    TemporalDelta Wt(t[tau], dt);

    SpatiallyConstant Ks(integrator.lengthInterval(),
                         integrator.spatialFreqResolution());
    TemporallyConstant Kt(integrator.timeInterval(),
                          integrator.temporalFreqResolution());


    cube F_e = zeros(r.n_elem, r.n_elem, t.n_elem);
    cx_cube G = zeros<cx_cube>(r.n_elem, r.n_elem, t.n_elem);

    for(int l=0; l < int(t.n_elem); l++){
        for(int i = 0; i < int(r.n_elem); i++){
            for(int j = 0; j < int(r.n_elem); j++){
                F_e(i,j,l) = 1.
                        / integrator.lengthInterval()
                        / integrator.lengthInterval()
                        /integrator.timeInterval();

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


    REQUIRE(diff_real.max() == Approx(0.0).epsilon(1e-7));
    REQUIRE(diff_imag.max() == Approx(0.0).epsilon(1e-7));

}


TEST_CASE("gaussConstConvolutionTest_test_0") {
    runGaussConstConvolutionTest(3, 0.05, 4, 0.05,
                                 2, 0.25);
}

TEST_CASE("gaussConstConvolutionTest_test_1") {
    runGaussConstConvolutionTest(2, 0.05, 3, 0.05,
                                 1,  0.25);
}
