#include "test_system_grc_irf_2.h"

test_system_grc_irf_2::test_system_grc_irf_2(string testLabel,
                                               string filename,
                                               double preCalls,
                                               double calls,
                                               double epsilon)
    : MCintegrationTest(testLabel, filename, preCalls, calls, epsilon)
{
    m_ndim=3;
}


void test_system_grc_irf_2::runTest()
{
    int ns = 9;
    int nt = 8;
    double dt = 2;
    double ds = 0.1;

    double w_g = 1.0;
    double w_rg = 1.0;
    double w_cr = 1.0;

    double a_rg = 0.1;

    double phaseDuration = 42.5 ;
    double dampingFactor = 0.38;

    double tau_rg=10;
    double tau_rc=39;

    double delay_rg=4;
    double delay_rc=20;


    double a_ri = 0.2;
    double a_rc= 0.1;

    double w_rc = 0.5;



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

    for(int i=0; i < nt-1; i++){
        m_t = integrator.timeVec()[int(pow(2,i))];

        //relayCell cell
        RelayCell relay(&integrator);

        SpatialGaussian Krg_s(a_rg);
        DecayingExponential Krg_t(tau_rg, delay_rg);
        SeparableKernel Krg(w_rg, &Krg_s, &Krg_t);

        SpatialGaussian Krc_s(a_rc);
        DecayingExponential Krc_t(tau_rc, delay_rc);
        SeparableKernel Krc(w_rc, &Krc_s, &Krc_t);


        //cortical cell
        CorticalCell cortical(&integrator);

        SpatialDelta Kcr_s(ds,{0,0});
        TemporalDelta Kcr_t(0, dt);
        SeparableKernel Kcr(w_cr, &Kcr_s, &Kcr_t);


        //Connect
        relay.addGanglionCell(&ganglion, Krg);
        relay.addCorticalCell(&cortical, Krc);
        cortical.addRelayCell(&relay, Kcr);


        if(m_computeMC){
            struct Param params = {this, stim, W, Krg, Krc, Kcr, m_peak};
            m_results.push_back(computeIntegral(xl, xu, &params));
        }

        relay.computeImpulseResponse();
        double ftIntegrator = relay.impulseResponse()
                (integrator.nPointsSpatial()/2,
                 integrator.nPointsSpatial()/2,
                 int(pow(2,i)));



        INFO(     "a_ri=" <<  a_ri
                  << "  a_rc=" << a_rc
                  << "  w_rc=" << w_rc
                  << "  t=" << m_t);
        CHECK(ftIntegrator == Approx(m_results[i]).epsilon(m_epsilon));
    }
    writeOutputFile();
}

double test_system_grc_irf_2::integrand(double *k, size_t dim, void *params)
{
    (void)(dim);
    double kx = k[0];
    double ky = k[1];
    double w = k[2];
    Param r = *(Param *) params;

    cx_double Wr = r.Wg.fourierTransform({kx, ky}, w)
            *(r.Krg.fourierTransform({kx, ky}, w))
            / (1. - r.Krc.fourierTransform({kx, ky}, w)
               * r.Kcr.fourierTransform({kx, ky}, w))
            * exp(-core::i*m_t*w);

    cx_double res =  Wr /(8*core::pi*core::pi*core::pi);
    return real(res);
}
