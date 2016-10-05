#include "test_system_gric_pg_1.h"

test_system_gric_pg_1::test_system_gric_pg_1(string testLabel,
                                             string filename,
                                             double preCalls,
                                             double calls,
                                             double epsilon)
    : MCintegrationTest(testLabel, filename, preCalls, calls, epsilon)
{
    m_ndim=2;
}


void test_system_gric_pg_1::runTest()
{
    int ns = 9;
    int nt = 1;
    double dt = 1;
    double ds = 0.1;

    double w_w = 1.0;
    double w_ig = 1.0;
    double w_rg = 1.0;
    double w_ri = -0.3;
    double w_cr = 1.0;
    vec weights = {0.0, 0.8};

    double C = 1.0;
    double orientation = 0.;
    int wId = 0;
    double phase = 0.0;
    vec kId = linspace(0, 80 , 5);
    vec maskSize = linspace(0.01, 13, 5);


    //integrator
    Integrator integrator(nt, dt, ns, ds);
    vec k = integrator.spatialFreqVec();
    vec w = integrator.temporalFreqVec();

    m_peak = integrator.temporalFreqResolution();
    double wd = w(wId);
    double xl[2] = {integrator.spatialFreqVec().min(),
                    integrator.spatialFreqVec().min()};
    double xu[2] ={integrator.spatialFreqVec().max(),
                   integrator.spatialFreqVec().max()};



    //ganglion cell
    DOG Ws(0.62, 1.26, 0.85);
    TemporalDelta Wt(0, dt);
    SeparableKernel W(w_w, &Ws, &Wt);
    GanglionCell ganglion(&integrator, W);


    int q = 0;
    for(double wc : weights){
        //relayCell cell
        RelayCell relay(&integrator);

        SpatialGaussian Krg_s(0.1);
        TemporalDelta Krg_t(0, dt);
        SeparableKernel Krg(w_rg, &Krg_s, &Krg_t);

        SpatialGaussian Kri_s(0.5);
        TemporalDelta Kri_t(0, dt);
        SeparableKernel Kri(w_ri, &Kri_s, &Kri_t);

        SpatialGaussian Krc_s(0.1);
        TemporalDelta Krc_t(0, dt);
        SeparableKernel Krc(wc, &Krc_s, &Krc_t);


        //interneuron cell
        Interneuron interneuron(&integrator);

        SpatialGaussian Kig_s(1.0);
        TemporalDelta Kig_t(0, dt);
        SeparableKernel Kig(w_ig, &Kig_s, &Kig_t);

        SpatialGaussian Kic_s(1.0);
        TemporalDelta Kic_t(0, dt);
        SeparableKernel Kic(wc, &Kic_s, &Kic_t);

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


        for(double d : maskSize){
            for(int i : kId){
                CircleMaskGrating stim(&integrator, k(i), wd,
                                       C, phase, orientation, d);
                stim.computeFourierTransform();

                if(m_computeMC){
                    struct Param params = {this, stim, W, Kig, Kic, Krg, Kri, Krc, Kcr, m_peak};
                    m_results.push_back(computeIntegral(xl, xu, &params));
                }

                relay.computeResponse(&stim, true);
                double ftIntegrator = relay.response()
                        (integrator.nPointsSpatial()/2,
                         integrator.nPointsSpatial()/2,
                         0);



                INFO( "kd=" <<  stim.spatialFreq()
                      << "  d=" << stim.maskSize()
                      << "   Ns=" << integrator.nPointsSpatial());
                CHECK(ftIntegrator == Approx(m_results[q]).epsilon(m_epsilon));
                q+=1;
            }

        }

    }
    writeOutputFile();
}

double test_system_gric_pg_1::integrand(double *k, size_t dim, void *params)
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

    cx_double res =  Wr *  r.S.fourierTransformAtFrequency({kx, ky}, wd)
            /(8*core::pi*core::pi*core::pi)*r.peak;
    return real(res);
}
