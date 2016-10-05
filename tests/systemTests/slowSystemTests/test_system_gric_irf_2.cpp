#include "test_system_gric_irf_2.h"

test_system_gric_irf_2::test_system_gric_irf_2(string testLabel,
                                               string filename,
                                               double preCalls,
                                               double calls,
                                               double epsilon)
    : MCintegrationTest(testLabel, filename, preCalls, calls, epsilon)
{
    m_ndim=3;
}


void test_system_gric_irf_2::runTest()
{
    int ns = 9;
    int nt = 8;
    double dt = 2;
    double ds = 0.1;

    double w_g = 1.0;
    double w_rg = 1.0;
    double w_cr = 1.0;

    double a_rg = 0.1;
    double a_ig = 0.3;

    double phaseDuration = 42.5 ;
    double dampingFactor = 0.38;

    double tau_rg=10;
    double tau_ri=10;
    double tau_rc=10;
    double tau_ig=10;
    double tau_ic=10;

    double delay_rg=0;
    double delay_ri=0;
    double delay_rc=0;
    double delay_ig=0;
    double delay_ic=0;



    double a_ri = 0.2;
    double a_rc= 0.1;
    double a_ic= 0.9;

    double w_ig = 1;
    double w_ri = -0.5;
    double w_rc = 0.5;
    double w_ic = 2.0;


    //integrator
    Integrator integrator(nt, dt, ns, ds);

    m_peak = integrator.temporalFreqResolution();
    double xl[3] = {integrator.spatialFreqVec().min(),
                    integrator.spatialFreqVec().min(),
                    integrator.temporalFreqVec().min()};
    double xu[3] ={integrator.spatialFreqVec().max(),
                   integrator.spatialFreqVec().max(),
                   integrator.temporalFreqVec().max()};


    //ganglion cell
    DOG Ws(0.62, 1.26, 0.85);
    Biphasic Wt(phaseDuration, dampingFactor, 0.0);
    SeparableKernel W(w_g, &Ws, &Wt);
    GanglionCell ganglion(&integrator, W);

    //stimulus
    CircleMaskGrating stim(&integrator,0, 0, 0, 0, 0, 0);

    int q = 0;


    //relayCell cell
    RelayCell relay(&integrator);

    SpatialGaussian Krg_s(a_rg);
    DecayingExponential Krg_t(tau_rg, delay_rg);
    SeparableKernel Krg(w_rg, &Krg_s, &Krg_t);

    SpatialGaussian Kri_s(a_ri);
    DecayingExponential Kri_t(tau_ri, delay_ri);
    SeparableKernel Kri(w_ri, &Kri_s, &Kri_t);

    SpatialGaussian Krc_s(a_rc);
    DecayingExponential Krc_t(tau_rc, delay_rc);
    SeparableKernel Krc(w_rc, &Krc_s, &Krc_t);


    //interneuron cell
    Interneuron interneuron(&integrator);

    SpatialGaussian Kig_s(a_ig);
    DecayingExponential Kig_t(tau_ig, delay_ig);
    SeparableKernel Kig(w_ig, &Kig_s, &Kig_t);

    SpatialGaussian Kic_s(a_ic);
    DecayingExponential Kic_t(tau_ic, delay_ic);
    SeparableKernel Kic(w_ic, &Kic_s, &Kic_t);

    //cortical cell
    CorticalCell cortical(&integrator);

    SpatialDelta Kcr_s(ds,{0,0});
    TemporalDelta Kcr_t(0, dt);
    SeparableKernel Kcr(w_cr, &Kcr_s, &Kcr_t);


    //Connect
    relay.addGanglionCell(&ganglion, Krg);
    relay.addInterNeuron(&interneuron, Kri);
    relay.addCorticalCell(&cortical, Krc);
    interneuron.addGanglionCell(&ganglion, Kig);
    interneuron.addCorticalCell(&cortical, Kic);
    cortical.addRelayCell(&relay, Kcr);


    if(m_computeMC){
        struct Param params = {this, stim, W, Kig, Kic, Krg, Kri, Krc, Kcr, m_peak};
        m_results.push_back(computeIntegral(xl, xu, &params));
    }

    relay.computeImpulseResponse();
    double ftIntegrator = relay.impulseResponse()
            (integrator.nPointsSpatial()/2,
             integrator.nPointsSpatial()/2,
             0);



    INFO(     "a_ri=" <<  a_ri
              << "  a_rc=" << a_rc
              << "  a_ic=" << a_ic
              << "  w_ig=" << w_ig
              << "  w_ri=" << w_ri
              << "  w_rc=" << w_rc
              << "  w_ic=" << w_ic);
    CHECK(ftIntegrator == Approx(m_results[q]).epsilon(m_epsilon));
    q+=1;
    writeOutputFile();
}

double test_system_gric_irf_2::integrand(double *k, size_t dim, void *params)
{
    (void)(dim);
    double kx = k[0];
    double ky = k[1];
    double w = k[2];
    Param r = *(Param *) params;

    cx_double Wr = r.Wg.fourierTransform({kx, ky}, w)
            *(r.Krg.fourierTransform({kx, ky}, w)
              + r.Kri.fourierTransform({kx, ky}, w)
              * r.Kig.fourierTransform({kx, ky}, w))
            / (1. -r. Kic.fourierTransform({kx, ky}, w)
               * r.Kcr.fourierTransform({kx, ky}, w)
               * r.Kri.fourierTransform({kx, ky}, w)
               - r.Krc.fourierTransform({kx, ky}, w)
               * r.Kcr.fourierTransform({kx, ky}, w));

    cx_double res =  Wr /(8*core::pi*core::pi*core::pi);
    return real(res);
}
