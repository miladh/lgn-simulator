/**********************************************************************
 *  Test: response of ganglion cell, relay cell and cortical cell with
 *  full-field grating stimulus:
 *
 *               Rg(r,t) = W(r,t) * S(r,t)
 *               Rc(r,t) = [Wr(r,t)Kcr(r,t)] * S(r,t)
 *               Rr(r,t) = [W(r,t)Krg(r,t)/ [1 - Krc(r,t)Kcr(r,t)] ] * S(r,t)
 *
 *  Analytic source: closed-form experssion
 *
 * ********************************************************************/

#include <unittest++/UnitTest++.h>
#include <lgnSimulator.h>

using namespace lgnSimulator;

void runSystemTest_GRC(int nt, double dt, int ns, double ds,
                       double C, int wdId, int kxId, int thetaId,
                       const Kernel &W, const Kernel &Krg,
                       const Kernel &Krc, const Kernel &Kcr)
{

    //integrator
    Integrator integrator(nt, dt, ns, ds);
    vec r = integrator.spatialVec();
    vec k = integrator.spatialFreqVec();
    vec t = integrator.timeVec();
    vec w = integrator.temporalFreqVec();
    vec orientations = vec{0., 90., 180., 270, -360.};


    //stimulus
    double wd = w(wdId);
    double spatialFreq = k(kxId);
    double orientation = orientations(thetaId);

    FullFieldGrating grating(integrator, spatialFreq, orientation, wd, C);
    grating.computeFourierTransform();

    //ganglion cell
    GanglionCell ganglion(integrator, W);

    //relayCell cell
    double R0 = 0.78;
    RelayCell relay(integrator, R0);

    //cortical cell
    CorticalCell cortical(integrator);

    //Connect
    relay.addGanglionCell(&ganglion, Krg);
    relay.addCorticalCell(&cortical, Krc);
    cortical.addRelayCell(&relay, Kcr);

    //Compute
    ganglion.computeResponse(&grating);
    relay.computeResponse(&grating);
    cortical.computeResponse(&grating);

    cube Rg_e = zeros(r.n_elem, r.n_elem, t.n_elem);
    cube Rr_e = zeros(r.n_elem, r.n_elem, t.n_elem);
    cube Rc_e = zeros(r.n_elem, r.n_elem, t.n_elem);

    double kx = spatialFreq*cos(orientation*core::pi/180.);
    double ky = spatialFreq*sin(orientation*core::pi/180.);
    complex<double> Wijl = W.fourierTransform({kx, ky}, wd);
    complex<double> Wr = Wijl* Krg.fourierTransform({kx, ky}, wd)
                        /(1. - Krc.fourierTransform({kx, ky}, wd)
                        *Kcr.fourierTransform({kx, ky}, wd));
    complex<double> Wc = Wr * Kcr.fourierTransform({kx, ky}, wd);

    for(int l=0; l < int(t.n_elem); l++){
        for(int i = 0; i < int(r.n_elem); i++){
            for(int j = 0; j < int(r.n_elem); j++){
                Rg_e(i,j,l) = C * abs(Wijl)
                        * cos(kx*r[i] + ky*r[j] - wd * t[l] + arg(Wijl));
                Rr_e(i,j,l) = R0 + C * abs(Wr)
                        * cos(kx*r[i] + ky*r[j] - wd * t[l] + arg(Wr));
                Rc_e(i,j,l) = C * abs(Wc)
                        * cos(kx*r[i] + ky*r[j] - wd * t[l] + arg(Wc));
            }
        }
    }

    cube diff_g = abs(Rg_e - ganglion.response());
    cube diff_r = abs(Rr_e - relay.response());
    cube diff_c = abs(Rc_e - cortical.response());

    CHECK_CLOSE(diff_g.max(), 0.0, 1e-10);
    CHECK_CLOSE(diff_r.max(), 0.0, 1e-10);
    CHECK_CLOSE(diff_c.max(), 0.0, 1e-10);


}



SUITE(system){

    TEST(runSystemTest_GRC_0){
        double dt = 0.1;
        double ds = 0.1;
        SpatialGaussian Ws(1, 0.25);
        TemporalDelta Wt(0, dt);
        SeparableKernel W(&Ws, &Wt);

        SpatialDelta Krg_s(0.5, ds, {0,0,});
        TemporalDelta Krg_t(2*dt, dt);
        SeparableKernel Krg(&Krg_s, &Krg_t);

        SpatialDelta Krc_s(0.5, ds, {0,0,});
        TemporalDelta Krc_t(0, dt);
        SeparableKernel Krc(&Krc_s, &Krc_t);

        SpatialGaussian Kcr_s(-1.5, 0.23);
        TemporalDelta Kcr_t(0, dt);
        SeparableKernel Kcr(&Kcr_s, &Kcr_t);

        runSystemTest_GRC(2, dt, 2, ds,
                          -1, 0, 1, 0,
                          W, Krg, Krc, Kcr);
    }


    TEST(runSystemTest_GRC_1){
        double dt = 0.1;
        double ds = 0.1;
        SpatialGaussian Ws(1, 0.25);
        TemporalDelta Wt(0, dt);
        SeparableKernel W(&Ws, &Wt);

        SpatialDelta Krg_s(0.5, ds, {0,0,});
        TemporalDelta Krg_t(2*dt, dt);
        SeparableKernel Krg(&Krg_s, &Krg_t);

        SpatialDelta Krc_s(0.5, ds, {0,0,});
        TemporalDelta Krc_t(0, dt);
        SeparableKernel Krc(&Krc_s, &Krc_t);

        SpatialGaussian Kcr_s(-1.0, 0.55);
        TemporalDelta Kcr_t(3*dt, dt);
        SeparableKernel Kcr(&Kcr_s, &Kcr_t);

        runSystemTest_GRC(3, dt, 4, 0.1,
                          -1, 3, 1, 4,
                         W, Krg, Krc, Kcr);
    }
}


