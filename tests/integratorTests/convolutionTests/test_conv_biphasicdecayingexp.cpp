/**********************************************************************
 *  Test: convolution theorem applied on F = W * K = ifft(fft(W)fft(K)),
 *        where W = Gauss(r) * decayingExp(t) and
 *        K = delta(r) * biphasic(t)
 *
 *  Analytic source: closed-form experssion (Sympy)
 *
 * ********************************************************************/

#include <lgnSimulator.h>
#include <catch.hpp>


using namespace lgnSimulator;


double definiteIntegral_db(double t, double limit,
                           double dampingFactor, double phaseDuration,
                           double tau_d, double delay) {

    double tau = limit;
    if(t - delay < 0.0){
        return 0.0;
    }else{
        if(t - delay - limit < 0){
            tau = t - delay;
        }

        double f_result = dampingFactor*phaseDuration*
                (-core::pi*tau_d*cos(core::pi*tau/phaseDuration)
                 + phaseDuration*sin(core::pi*tau/phaseDuration))
                 * exp((delay - t + tau)/tau_d)/(pow(core::pi, 2) *pow(tau_d, 2)
                 + pow(phaseDuration, 2));
        return f_result;
    }



}


double convBiphasicDecay(double t,
                         double dampingFactor, double phaseDuration,
                         double tau_d, double delay){

    double term_1_lim_0 = definiteIntegral_db(t, 0, 1, phaseDuration, tau_d, delay);
    double term_1_lim_a = definiteIntegral_db(t, phaseDuration, 1, phaseDuration, tau_d, delay);
    double term_2_lim_a = definiteIntegral_db(t, phaseDuration, dampingFactor, phaseDuration, tau_d, delay)
                        * Special::heaviside(t - delay - phaseDuration);
    double term_2_lim_2a = definiteIntegral_db(t, 2*phaseDuration,dampingFactor, phaseDuration, tau_d, delay)
                         * Special::heaviside(t - delay - phaseDuration);

    return (term_1_lim_a - term_1_lim_0) + (term_2_lim_2a - term_2_lim_a);
}


void runDogdecayingExpConvolutionTest(int nt, double dt, int ns, double ds,
                                      double phaseDuration, double dampingFactor,
                                      double a, double tau_d, int delay,
                                      vec2 shift)
{

    Integrator integrator(nt, dt, ns, ds);
    vec r = integrator.spatialVec();
    vec k = integrator.spatialFreqVec();
    vec t = integrator.timeVec();
    vec w = integrator.temporalFreqVec();

    SpatialGaussian Ws( a);
    Biphasic Wt(phaseDuration, dampingFactor, 0);

    SpatialDelta Ks( ds, shift);
    DecayingExponential Kt(tau_d, t[delay]);

    cube F_e = zeros(r.n_elem, r.n_elem, t.n_elem);
    cx_cube G = zeros<cx_cube>(r.n_elem, r.n_elem, t.n_elem);

    for(int l=0; l < int(t.n_elem); l++){
        for(int i = 0; i < int(r.n_elem); i++){
            for(int j = 0; j < int(r.n_elem); j++){
                F_e(i,j,l) =  Ws.spatial(vec2{r[i], r[j]} - shift)
                        * convBiphasicDecay(t[l], dampingFactor, phaseDuration,
                                            tau_d, t[delay]);

                G(i,j,l) = Ws.fourierTransform({k[i], k[j]})
                        * Wt.fourierTransform(w[l])
                        * Ks.fourierTransform({k[i], k[j]})
                        * Kt.fourierTransform(w[l]);
            }
        }
    }

    cube F = integrator.backwardFFT(G);
    cube diff = (F_e-F)*ds*ds; // divide by the contributions from spatial part

    cube diff_real = abs(real(diff));
    cube diff_imag = abs(imag(diff));

    //    cout << diff_real.max() << endl;
    //    cout << diff_imag.max() << endl;

    REQUIRE(diff_real.max() == Approx(0.0).epsilon(1e-5));
    REQUIRE(diff_imag.max() == Approx(0.0).epsilon(1e-5));

}



TEST_CASE("runDogdecayingExpConvolutionTest_test_0") {
    double phaseDuration = 42.5;
    double dampingFactor = 0.38;
    double a = 0.25;
    double tau_d = 1.0;
    int delay=0;

    runDogdecayingExpConvolutionTest(8, 1.0, 5, 0.05,
                                     phaseDuration,dampingFactor,
                                     a, tau_d,  delay,
                                     vec2{0.0, 0.0});
}


TEST_CASE("runDogdecayingExpConvolutionTest_test_1") {
    double phaseDuration = 42.5;
    double dampingFactor = 0.38;
    double a = 0.25;
    double tau_d = 22;
    int delay=10;

    runDogdecayingExpConvolutionTest(8, 1.0, 5, 0.05,
                                     phaseDuration,dampingFactor,
                                     a, tau_d,  delay,
                                     vec2{0.0, 0.0});
}

