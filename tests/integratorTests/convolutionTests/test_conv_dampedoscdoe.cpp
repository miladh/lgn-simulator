/**********************************************************************
 *  Test: convolution theorem applied on F = W * K = ifft(fft(W)fft(K)),
 *        where W = delta(r) * dampedOscillator(t) and
 *        K = delta(r) * combinedRC(t)
 *
 *  Analytic source: closed-form experssion
 *
 * ********************************************************************/

#include <unittest++/UnitTest++.h>
#include <lgnSimulator.h>

using namespace lgnSimulator;


double convDampedOscDOE(double phaseDuration, double dampingFactor, int delay_osc,
                               double cenLatency, double surLatency, int delay_rc,
                               double t)
{
    double f_result;
       f_result = pow(phaseDuration, 2)*(pow(core::pi, 4)*dampingFactor*delay_rc*pow(surLatency, 4)*sin(core::pi*t/phaseDuration) - 2*pow(core::pi, 4)*dampingFactor*phaseDuration*pow(surLatency, 4)*sin(core::pi*t/phaseDuration) + 2*pow(core::pi, 3)*dampingFactor*delay_rc*phaseDuration*pow(surLatency, 3)*cos(core::pi*t/phaseDuration) - 4*pow(core::pi, 3)*dampingFactor*pow(phaseDuration, 2)*pow(surLatency, 3)*cos(core::pi*t/phaseDuration) + 2*pow(core::pi, 3)*dampingFactor*phaseDuration*pow(surLatency, 4)*cos(core::pi*t/phaseDuration) - 6*pow(core::pi, 2)*dampingFactor*pow(phaseDuration, 2)*pow(surLatency, 3)*sin(core::pi*t/phaseDuration) + 2*core::pi*dampingFactor*delay_rc*pow(phaseDuration, 3)*surLatency*cos(core::pi*t/phaseDuration) - 4*core::pi*dampingFactor*pow(phaseDuration, 4)*surLatency*cos(core::pi*t/phaseDuration) - 6*core::pi*dampingFactor*pow(phaseDuration, 3)*pow(surLatency, 2)*cos(core::pi*t/phaseDuration) - dampingFactor*delay_rc*pow(phaseDuration, 4)*sin(core::pi*t/phaseDuration) + 2*dampingFactor*pow(phaseDuration, 5)*sin(core::pi*t/phaseDuration) + 2*dampingFactor*pow(phaseDuration, 4)*surLatency*sin(core::pi*t/phaseDuration) + (-pow(core::pi, 4)*delay_rc*pow(surLatency, 4)*sin(core::pi*t/phaseDuration) - 2*pow(core::pi, 3)*delay_rc*phaseDuration*pow(surLatency, 3)*cos(core::pi*t/phaseDuration) - 2*pow(core::pi, 3)*phaseDuration*pow(surLatency, 4)*cos(core::pi*t/phaseDuration) + 6*pow(core::pi, 2)*pow(phaseDuration, 2)*pow(surLatency, 3)*sin(core::pi*t/phaseDuration) - 2*core::pi*delay_rc*pow(phaseDuration, 3)*surLatency*cos(core::pi*t/phaseDuration) + 6*core::pi*pow(phaseDuration, 3)*pow(surLatency, 2)*cos(core::pi*t/phaseDuration) + delay_rc*pow(phaseDuration, 4)*sin(core::pi*t/phaseDuration) - 2*pow(phaseDuration, 4)*surLatency*sin(core::pi*t/phaseDuration))*exp(2*phaseDuration/surLatency) + (pow(core::pi, 4)*dampingFactor*delay_rc*pow(surLatency, 4)*sin(core::pi*t/phaseDuration) - pow(core::pi, 4)*dampingFactor*phaseDuration*pow(surLatency, 4)*sin(core::pi*t/phaseDuration) - pow(core::pi, 4)*delay_rc*pow(surLatency, 4)*sin(core::pi*t/phaseDuration) + pow(core::pi, 4)*phaseDuration*pow(surLatency, 4)*sin(core::pi*t/phaseDuration) + 2*pow(core::pi, 3)*dampingFactor*delay_rc*phaseDuration*pow(surLatency, 3)*cos(core::pi*t/phaseDuration) - 2*pow(core::pi, 3)*dampingFactor*pow(phaseDuration, 2)*pow(surLatency, 3)*cos(core::pi*t/phaseDuration) + 2*pow(core::pi, 3)*dampingFactor*phaseDuration*pow(surLatency, 4)*cos(core::pi*t/phaseDuration) - 2*pow(core::pi, 3)*delay_rc*phaseDuration*pow(surLatency, 3)*cos(core::pi*t/phaseDuration) + 2*pow(core::pi, 3)*pow(phaseDuration, 2)*pow(surLatency, 3)*cos(core::pi*t/phaseDuration) - 2*pow(core::pi, 3)*phaseDuration*pow(surLatency, 4)*cos(core::pi*t/phaseDuration) - 6*pow(core::pi, 2)*dampingFactor*pow(phaseDuration, 2)*pow(surLatency, 3)*sin(core::pi*t/phaseDuration) + 6*pow(core::pi, 2)*pow(phaseDuration, 2)*pow(surLatency, 3)*sin(core::pi*t/phaseDuration) + 2*core::pi*dampingFactor*delay_rc*pow(phaseDuration, 3)*surLatency*cos(core::pi*t/phaseDuration) - 2*core::pi*dampingFactor*pow(phaseDuration, 4)*surLatency*cos(core::pi*t/phaseDuration) - 6*core::pi*dampingFactor*pow(phaseDuration, 3)*pow(surLatency, 2)*cos(core::pi*t/phaseDuration) - 2*core::pi*delay_rc*pow(phaseDuration, 3)*surLatency*cos(core::pi*t/phaseDuration) + 2*core::pi*pow(phaseDuration, 4)*surLatency*cos(core::pi*t/phaseDuration) + 6*core::pi*pow(phaseDuration, 3)*pow(surLatency, 2)*cos(core::pi*t/phaseDuration) - dampingFactor*delay_rc*pow(phaseDuration, 4)*sin(core::pi*t/phaseDuration) + dampingFactor*pow(phaseDuration, 5)*sin(core::pi*t/phaseDuration) + 2*dampingFactor*pow(phaseDuration, 4)*surLatency*sin(core::pi*t/phaseDuration) + delay_rc*pow(phaseDuration, 4)*sin(core::pi*t/phaseDuration) - pow(phaseDuration, 5)*sin(core::pi*t/phaseDuration) - 2*pow(phaseDuration, 4)*surLatency*sin(core::pi*t/phaseDuration))*exp(phaseDuration/surLatency))*exp(delay_rc/surLatency)*exp(-2*phaseDuration/surLatency)/pow(pow(core::pi, 2)*pow(surLatency, 2) + pow(phaseDuration, 2), 3) - pow(phaseDuration, 2)*(pow(cenLatency, 4)*pow(core::pi, 4)*dampingFactor*delay_rc*sin(core::pi*t/phaseDuration) - 2*pow(cenLatency, 4)*pow(core::pi, 4)*dampingFactor*phaseDuration*sin(core::pi*t/phaseDuration) + 2*pow(cenLatency, 4)*pow(core::pi, 3)*dampingFactor*phaseDuration*cos(core::pi*t/phaseDuration) + 2*pow(cenLatency, 3)*pow(core::pi, 3)*dampingFactor*delay_rc*phaseDuration*cos(core::pi*t/phaseDuration) - 4*pow(cenLatency, 3)*pow(core::pi, 3)*dampingFactor*pow(phaseDuration, 2)*cos(core::pi*t/phaseDuration) - 6*pow(cenLatency, 3)*pow(core::pi, 2)*dampingFactor*pow(phaseDuration, 2)*sin(core::pi*t/phaseDuration) - 6*pow(cenLatency, 2)*core::pi*dampingFactor*pow(phaseDuration, 3)*cos(core::pi*t/phaseDuration) + 2*cenLatency*core::pi*dampingFactor*delay_rc*pow(phaseDuration, 3)*cos(core::pi*t/phaseDuration) - 4*cenLatency*core::pi*dampingFactor*pow(phaseDuration, 4)*cos(core::pi*t/phaseDuration) + 2*cenLatency*dampingFactor*pow(phaseDuration, 4)*sin(core::pi*t/phaseDuration) - dampingFactor*delay_rc*pow(phaseDuration, 4)*sin(core::pi*t/phaseDuration) + 2*dampingFactor*pow(phaseDuration, 5)*sin(core::pi*t/phaseDuration) + (-pow(cenLatency, 4)*pow(core::pi, 4)*delay_rc*sin(core::pi*t/phaseDuration) - 2*pow(cenLatency, 4)*pow(core::pi, 3)*phaseDuration*cos(core::pi*t/phaseDuration) - 2*pow(cenLatency, 3)*pow(core::pi, 3)*delay_rc*phaseDuration*cos(core::pi*t/phaseDuration) + 6*pow(cenLatency, 3)*pow(core::pi, 2)*pow(phaseDuration, 2)*sin(core::pi*t/phaseDuration) + 6*pow(cenLatency, 2)*core::pi*pow(phaseDuration, 3)*cos(core::pi*t/phaseDuration) - 2*cenLatency*core::pi*delay_rc*pow(phaseDuration, 3)*cos(core::pi*t/phaseDuration) - 2*cenLatency*pow(phaseDuration, 4)*sin(core::pi*t/phaseDuration) + delay_rc*pow(phaseDuration, 4)*sin(core::pi*t/phaseDuration))*exp(2*phaseDuration/cenLatency) + (pow(cenLatency, 4)*pow(core::pi, 4)*dampingFactor*delay_rc*sin(core::pi*t/phaseDuration) - pow(cenLatency, 4)*pow(core::pi, 4)*dampingFactor*phaseDuration*sin(core::pi*t/phaseDuration) - pow(cenLatency, 4)*pow(core::pi, 4)*delay_rc*sin(core::pi*t/phaseDuration) + pow(cenLatency, 4)*pow(core::pi, 4)*phaseDuration*sin(core::pi*t/phaseDuration) + 2*pow(cenLatency, 4)*pow(core::pi, 3)*dampingFactor*phaseDuration*cos(core::pi*t/phaseDuration) - 2*pow(cenLatency, 4)*pow(core::pi, 3)*phaseDuration*cos(core::pi*t/phaseDuration) + 2*pow(cenLatency, 3)*pow(core::pi, 3)*dampingFactor*delay_rc*phaseDuration*cos(core::pi*t/phaseDuration) - 2*pow(cenLatency, 3)*pow(core::pi, 3)*dampingFactor*pow(phaseDuration, 2)*cos(core::pi*t/phaseDuration) - 2*pow(cenLatency, 3)*pow(core::pi, 3)*delay_rc*phaseDuration*cos(core::pi*t/phaseDuration) + 2*pow(cenLatency, 3)*pow(core::pi, 3)*pow(phaseDuration, 2)*cos(core::pi*t/phaseDuration) - 6*pow(cenLatency, 3)*pow(core::pi, 2)*dampingFactor*pow(phaseDuration, 2)*sin(core::pi*t/phaseDuration) + 6*pow(cenLatency, 3)*pow(core::pi, 2)*pow(phaseDuration, 2)*sin(core::pi*t/phaseDuration) - 6*pow(cenLatency, 2)*core::pi*dampingFactor*pow(phaseDuration, 3)*cos(core::pi*t/phaseDuration) + 6*pow(cenLatency, 2)*core::pi*pow(phaseDuration, 3)*cos(core::pi*t/phaseDuration) + 2*cenLatency*core::pi*dampingFactor*delay_rc*pow(phaseDuration, 3)*cos(core::pi*t/phaseDuration) - 2*cenLatency*core::pi*dampingFactor*pow(phaseDuration, 4)*cos(core::pi*t/phaseDuration) - 2*cenLatency*core::pi*delay_rc*pow(phaseDuration, 3)*cos(core::pi*t/phaseDuration) + 2*cenLatency*core::pi*pow(phaseDuration, 4)*cos(core::pi*t/phaseDuration) + 2*cenLatency*dampingFactor*pow(phaseDuration, 4)*sin(core::pi*t/phaseDuration) - 2*cenLatency*pow(phaseDuration, 4)*sin(core::pi*t/phaseDuration) - dampingFactor*delay_rc*pow(phaseDuration, 4)*sin(core::pi*t/phaseDuration) + dampingFactor*pow(phaseDuration, 5)*sin(core::pi*t/phaseDuration) + delay_rc*pow(phaseDuration, 4)*sin(core::pi*t/phaseDuration) - pow(phaseDuration, 5)*sin(core::pi*t/phaseDuration))*exp(phaseDuration/cenLatency))*exp(delay_rc/cenLatency)*exp(-2*phaseDuration/cenLatency)/pow(pow(cenLatency, 2)*pow(core::pi, 2) + pow(phaseDuration, 2), 3);
       return f_result;
}

void runDampedOscCombinedRCConvolutionTest(int nt, double dt, int ns, double ds,
                                           double phaseDuration,
                                           double dampingFactor,
                                           int delay_osc,
                                           double wg, vec2 rg,
                                           double cenLatency,
                                           double surLatency,
                                           int delay_rc,
                                           double wk, vec2 rk)
{

    Integrator integrator(nt, dt, ns, ds);
    vec r = integrator.spatialVec();
    vec k = integrator.spatialFreqVec();
    vec t = integrator.timeVec();
    vec w = integrator.temporalFreqVec();

    SpatialDelta Ws(wg, ds, rg);
    DampedOscillator Wt(phaseDuration, dampingFactor, t[delay_osc]);

    SpatialDelta Ks(wk, ds, rk);
    DOE Kt(cenLatency, surLatency, t[delay_rc]);

    cube F_e = zeros(r.n_elem, r.n_elem, t.n_elem);
    cx_cube G = zeros<cx_cube>(r.n_elem, r.n_elem, t.n_elem);

    for(int l=0; l < int(t.n_elem); l++){
        for(int i = 0; i < int(r.n_elem); i++){
            for(int j = 0; j < int(r.n_elem); j++){
                F_e(i,j,l) = wk
                        * Ws.spatial(vec2{r[i], r[j]} - rk)
                        * convDampedOscDOE(phaseDuration, dampingFactor,0,
                                                  cenLatency, surLatency, t[delay_rc],
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

//    cout << diff_real << endl;
//    cout << real(F)/wk/wg*ds*ds << endl;

    CHECK_CLOSE(diff_real.max(), 0.0, 1e-3);
    CHECK_CLOSE(diff_imag.max(), 0.0, 1e-12);

}


SUITE(integrator){


//    TEST(dampedOscDOEConvolutionTest_test_0) {
//        runDampedOscCombinedRCConvolutionTest(8, 0.001, 2, 0.05,
//                                              42.5, 0.38, 0,
//                                              1., vec2{0.0, 0.0},
//                                              16., 32., 0,
//                                              1.0, vec2{0.0, 0.0});
//    }

}
