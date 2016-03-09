/**********************************************************************
 *  Test: convolution theorem applied on F = W * K = ifft(fft(W)fft(K)),
 *        where W = delta(r) * dampedOscillator(t) and
 *        K = delta(r) * doe(t)
 *
 *  Analytic source: closed-form experssion
 *
 * ********************************************************************/

#include <unittest++/UnitTest++.h>
#include <lgnSimulator.h>

using namespace lgnSimulator;


double definiteIntegral(double t, double limit, double A, double a,
                        double gamma, double delay)
{
    double tau = limit;
    if(t - delay < 0.0){
        return 0.0;
    }else{
        if(t - delay - limit < 0){
            tau = t - delay;
        }

        double factor = A/(a*a + core::pi*core::pi * gamma*gamma)
                /(a*a + core::pi*core::pi * gamma*gamma);
        double expTerm = exp((delay-t+tau)/gamma);
        double trigFactor = (core::pi * tau)/a;
        double cosTerm = core::pi*a*cos(trigFactor)
                *(a*a * (-2.*gamma + delay - t + tau)
                  + core::pi*core::pi * gamma*gamma * (delay-t+tau) );
        double sinTerm = a*a/gamma * sin(trigFactor)
                *(a*a *(-gamma + delay - t + tau)
                  + core::pi*core::pi * gamma*gamma *(gamma + delay - t + tau));

        return factor * expTerm * (cosTerm - sinTerm);
    }
}


double convDampedOscDOE(double phaseDuration, double dampingFactor,
                        double cenLatency, double surLatency, double delay,
                        double t){

    double cenTerm1_1 = definiteIntegral(t, 0, 1, phaseDuration, cenLatency, delay);
    double cenTerm1_2 = definiteIntegral(t, phaseDuration, 1, phaseDuration, cenLatency, delay);
    double cenTerm2_1 = definiteIntegral(t, phaseDuration, dampingFactor, phaseDuration, cenLatency,
                                         delay)*Special::heaviside(t - delay - phaseDuration);
    double cenTerm2_2 = definiteIntegral(t, 2*phaseDuration,dampingFactor, phaseDuration, cenLatency,
                                         delay)*Special::heaviside(t - delay - phaseDuration);

    double cenTerm = (cenTerm1_2 - cenTerm1_1) + (cenTerm2_2 - cenTerm2_1);

    double surTerm1_1 = definiteIntegral(t, 0, 1., phaseDuration,  surLatency, delay);
    double surTerm1_2 = definiteIntegral(t, phaseDuration, 1., phaseDuration, surLatency, delay);
    double surTerm2_1 = definiteIntegral(t, phaseDuration, dampingFactor, phaseDuration, surLatency,
                                         delay)*Special::heaviside(t - delay - phaseDuration);
    double surTerm2_2 = definiteIntegral(t, 2*phaseDuration, dampingFactor, phaseDuration,  surLatency,
                                         delay)*Special::heaviside(t - delay - phaseDuration);


    double surTerm = (surTerm1_2 - surTerm1_1) + (surTerm2_2 - surTerm2_1);

    return cenTerm - surTerm;
}

void runDampedOscCombinedRCConvolutionTest(int nt, double dt, int ns, double ds,
                                           double phaseDuration,
                                           double dampingFactor,
                                           double wg, vec2 rg,
                                           double cenLatency,
                                           double surLatency,
                                           int delay_doe,
                                           double wk, vec2 rk)
{

    Integrator integrator(nt, dt, ns, ds);
    vec r = integrator.spatialVec();
    vec k = integrator.spatialFreqVec();
    vec t = integrator.timeVec();
    vec w = integrator.temporalFreqVec();

    SpatialDelta Ws(wg, ds, rg);
    Biphasic Wt(phaseDuration, dampingFactor, 0);

    SpatialDelta Ks(wk, ds, rk);
    DOE Kt(cenLatency, surLatency, t[delay_doe]);

    cube F_e = zeros(r.n_elem, r.n_elem, t.n_elem);
    cx_cube G = zeros<cx_cube>(r.n_elem, r.n_elem, t.n_elem);

    for(int l=0; l < int(t.n_elem); l++){
        for(int i = 0; i < int(r.n_elem); i++){
            for(int j = 0; j < int(r.n_elem); j++){
                F_e(i,j,l) = wk
                        * Ws.spatial(vec2{r[i], r[j]} - rk)
                        * convDampedOscDOE(phaseDuration, dampingFactor,
                                           cenLatency, surLatency, t[delay_doe],
                                           t[l]);

                G(i,j,l) = Ws.fourierTransform({k[i], k[j]})
                        * Wt.fourierTransform(w[l])
                        * Ks.fourierTransform({k[i], k[j]})
                        * Kt.fourierTransform(w[l]);
            }
        }
    }

    cx_cube F = integrator.backwardFFT(G);
    cx_cube diff = (F_e-F)/wk/wg*ds*ds; // divide by the contributions from spatial part

    cube diff_real = abs(real(diff));
    cube diff_imag = abs(imag(diff));

//    cout << diff_real.max() << endl;
//    cout << diff_imag.max() << endl;


    CHECK_CLOSE(diff_real.max(), 0.0, 1e-3);
    CHECK_CLOSE(diff_imag.max(), 0.0, 1e-12);

}


SUITE(integrator){


        TEST(dampedOscDOEConvolutionTest_test_0) {
            runDampedOscCombinedRCConvolutionTest(12, 0.1, 2, 0.05,
                                                  42.5, 0.38,
                                                  1., vec2{0.0, 0.0},
                                                  16., 32., 0,
                                                  1.0, vec2{0.0, 0.0});
        }



        TEST(dampedOscDOEConvolutionTest_test_1) {
            runDampedOscCombinedRCConvolutionTest(12, 0.1, 2, 0.05,
                                                  42.5, 0.38,
                                                  1., vec2{0.0, 0.0},
                                                  16., 32., 1365,
                                                  1.0, vec2{0.0, 0.0});
        }

}
