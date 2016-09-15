#include "test_system_gric_irf_1.h"

test_system_gric_irf_1::test_system_gric_irf_1(string testLabel,
                                               string filename,
                                               double preCalls,
                                               double calls,
                                               double epsilon)
    : MCintegrationTest(testLabel, filename, preCalls, calls, epsilon)
{

}


void test_system_gric_irf_1::runTest()
{
    int ns = 9;
    int nt = 1;
    double dt = 1;
    double ds = 0.1;

    double w_g = 1.0;
    double w_rg = 1.0;
    double w_cr = 1.0;

    double a_rg = 0.1;
    double a_ig = 0.3;


    vec a_ri_vec = linspace(0.1, 2.5, 2);
    vec a_rc_vec= linspace(0.1, 2.5, 2);
    vec a_ic_vec= linspace(0.1, 2.5, 2);

    vec w_ig_vec = linspace(0, 1, 2);
    vec w_ri_vec = linspace(-4, 0, 2);
    vec w_rc_vec = linspace(0, 0.9, 2);
    vec w_ic_vec = linspace(0, 4, 2);


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
    for(double a_ri : a_ri_vec){
        for(double a_rc : a_rc_vec){
            for(double a_ic : a_ic_vec){
                for(double w_ig : w_ig_vec){
                    for(double w_ri : w_ri_vec){
                        for(double w_rc : w_rc_vec){
                            for(double w_ic : w_ic_vec){

                                //relayCell cell
                                RelayCell relay(&integrator);

                                SpatialGaussian Krg_s(a_rg);
                                TemporalDelta Krg_t(0, dt);
                                SeparableKernel Krg(w_rg, &Krg_s, &Krg_t);

                                SpatialGaussian Kri_s(a_ri);
                                TemporalDelta Kri_t(0, dt);
                                SeparableKernel Kri(w_ri, &Kri_s, &Kri_t);

                                SpatialGaussian Krc_s(a_rc);
                                TemporalDelta Krc_t(0, dt);
                                SeparableKernel Krc(w_rc, &Krc_s, &Krc_t);


                                //interneuron cell
                                Interneuron interneuron(&integrator);

                                SpatialGaussian Kig_s(a_ig);
                                TemporalDelta Kig_t(0, dt);
                                SeparableKernel Kig(w_ig, &Kig_s, &Kig_t);

                                SpatialGaussian Kic_s(a_ic);
                                TemporalDelta Kic_t(0, dt);
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
                            }

                        }
                    }
                }
            }
        }
    }
    writeOutputFile();
}

double test_system_gric_irf_1::integrand(double *k, size_t dim, void *params)
{
    (void)(dim);
    double wd = 0.0;
    double kx = k[0];
    double ky = k[1];
    Param r = *(Param *) params;

    cx_double Wr = r.Wg.fourierTransform({kx, ky}, wd)
            *(r.Krg.fourierTransform({kx, ky}, wd)
              + r.Kri.fourierTransform({kx, ky}, wd)
              * r.Kig.fourierTransform({kx, ky}, wd))
            / (1. -r. Kic.fourierTransform({kx, ky}, wd)
               * r.Kcr.fourierTransform({kx, ky}, wd)
               * r.Kri.fourierTransform({kx, ky}, wd)
               - r.Krc.fourierTransform({kx, ky}, wd)
               * r.Kcr.fourierTransform({kx, ky}, wd));

    cx_double res =  Wr /(4*core::pi*core::pi)/**r.peak*/;
    return real(res);
}
