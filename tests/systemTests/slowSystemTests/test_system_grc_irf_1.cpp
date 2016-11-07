#include "test_system_grc_irf_1.h"

test_system_grc_irf_1::test_system_grc_irf_1(string testLabel,
                                               string filename,
                                               double preCalls,
                                               double calls,
                                               double epsilon)
    : MCintegrationTest(testLabel, filename, preCalls, calls, epsilon)
{
    m_ndim=2;
}


void test_system_grc_irf_1::runTest()
{
    int ns = 9;
    int nt = 1;
    double dt = 1;
    double ds = 0.1;

    double w_g = 1.0;
    double w_rg = 1.0;
    double w_cr = 1.0;

    double a_rg = 0.1;


    vec a_rc_vec= linspace(0.1, 2.5, 2);
    vec w_rc_vec = linspace(0, 0.9, 2);


    //integrator
    Integrator integrator(nt, dt, ns, ds);

    m_peak = integrator.temporalFreqResolution();
    double xl[2] = {integrator.spatialFreqVec().min(),
                    integrator.spatialFreqVec().min()};
    double xu[2] ={integrator.spatialFreqVec().max(),
                   integrator.spatialFreqVec().max()};



    //ganglion cell
    DOG Ws(0.62, 1.26, 0.85);
    TemporalDelta Wt(0, dt);
    SeparableKernel W(w_g, &Ws, &Wt);
    GanglionCell ganglion(&integrator, W);


    //stimulus
    CircleMaskGrating stim(&integrator,0, 0, 0, 0, 0, 0);

    int q = 0;
        for(double a_rc : a_rc_vec){
            for(double w_rc : w_rc_vec){

                //relayCell cell
                RelayCell relay(&integrator);

                SpatialGaussian Krg_s(a_rg);
                TemporalDelta Krg_t(0, dt);
                SeparableKernel Krg(w_rg, &Krg_s, &Krg_t);


                SpatialGaussian Krc_s(a_rc);
                TemporalDelta Krc_t(0, dt);
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
                         0);



                INFO(  "  a_rc=" << a_rc
                      << "  w_rc=" << w_rc);
                CHECK(ftIntegrator == Approx(m_results[q]).epsilon(m_epsilon));
                q+=1;
            }
    }
    writeOutputFile();
}

double test_system_grc_irf_1::integrand(double *k, size_t dim, void *params)
{
    (void)(dim);
    double wd = 0.0;
    double kx = k[0];
    double ky = k[1];
    Param r = *(Param *) params;

    cx_double Wr = r.Wg.fourierTransform({kx, ky}, wd)
            *(r.Krg.fourierTransform({kx, ky}, wd))
            / (1. - r.Krc.fourierTransform({kx, ky}, wd)
               * r.Kcr.fourierTransform({kx, ky}, wd));

    cx_double res =  Wr /(4*core::pi*core::pi)/**r.peak*/;
    return real(res);
}
