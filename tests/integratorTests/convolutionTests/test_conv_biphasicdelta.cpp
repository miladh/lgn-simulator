/**********************************************************************
 *  Test: convolution theorem applied on F = W * K = ifft(fft(W)fft(K)),
 *        where W = delta(r) * dampedOscillator(t) and
 *        K = delta(r) * delta(t)
 *
 *  Analytic source: closed-form experssion
 *
 * ********************************************************************/

#include <unittest++/UnitTest++.h>
#include <lgnSimulator.h>

using namespace lgnSimulator;


void runDampedOscConvolutionTest(int nt, double dt, int ns, double ds,
                                 double phaseDuration, double dampedFactor,
                                 int delay_osc,
                                 vec2 rg,
                                 int delay_id,  vec2 rk)
{

    Integrator integrator(nt, dt, ns, ds);
    vec r = integrator.spatialVec();
    vec k = integrator.spatialFreqVec();
    vec t = integrator.timeVec();
    vec w = integrator.temporalFreqVec();

    SpatialDelta Ws(ds, rg);
    Biphasic Wt(phaseDuration, dampedFactor, t[delay_osc]);

    SpatialDelta Ks( ds, rk);
    TemporalDelta Kt(t[delay_id], dt);

    cube F_e = zeros(r.n_elem, r.n_elem, t.n_elem);
    cx_cube G = zeros<cx_cube>(r.n_elem, r.n_elem, t.n_elem);

    for(int l=0; l < int(t.n_elem); l++){
        for(int i = 0; i < int(r.n_elem); i++){
            for(int j = 0; j < int(r.n_elem); j++){
                F_e(i,j,l) =  Ws.spatial(vec2{r[i], r[j]} - rk)
                        * Wt.temporal(t[l] - t[delay_id]);

                G(i,j,l) = Ws.fourierTransform({k[i], k[j]})
                        * Wt.fourierTransform(w[l])
                        * Ks.fourierTransform({k[i], k[j]})
                        * Kt.fourierTransform(w[l]);
            }
        }
    }

    cx_cube F = integrator.backwardFFT(G);
    cx_cube diff = (F - F_e)*ds*ds; // divide by the contributions from spatial part

    cube diff_real = abs(real(diff));
    cube diff_imag = abs(imag(diff));

//    cout << diff_real.max() << endl;
//    cout << diff_imag.max() << endl;

    CHECK_CLOSE(diff_real.max(), 0.0, 1e-3);
    CHECK_CLOSE(diff_imag.max(), 0.0, 1e-12);

}


SUITE(integrator){


    TEST(runDampedOscConvolutionTest_test_0) {
        runDampedOscConvolutionTest(10, 0.1, 2, 0.05,
                                    42.5, 0.38, 0.0,
                                    vec2{0.0, 0.0},
                                    0,  vec2{0.0, 0.0});
    }

    TEST(runDampedOscConvolutionTest_test_1) {
        runDampedOscConvolutionTest(10, 0.1, 2, 0.05,
                                    42.5, 0.38, 0.0,
                                    vec2{0.0, 0.0},
                                    32, vec2{0.0, 0.0});
    }

    TEST(runDampedOscConvolutionTest_test_2) {
        runDampedOscConvolutionTest(10, 0.1, 2, 0.05,
                                    42.5, 0.38, 20,
                                    vec2{0.0, 0.0},
                                    50, vec2{0.05, 0.05});
    }


}
