/**********************************************************************
 *  Test: response of ganglion cell, relay cell, interneuron and
 *  cortical cell with full-field grating stimulus:
 *
 *               Rg(r,t) = Wg(r,t) * S(r,t)
 *               Ri(r,t) = [Wg(r,t)Kig(r,t) + Kic(r,t)Kcr(r,t)*Wr(r,t)] * S(r,t)
 *               Rc(r,t) = [Wr(r,t)Kcr(r,t)] * S(r,t)
 *               Rr(r,t) = [[W(r,t)Krg(r,t) + Kri(r,t) Wg(r,t)Kig(r,t)]
 *                         /[1 - Kic(r,t)Kcr(r,t)] - Krc(r,t)Kcr(r,t)] ] * S(r,t)
 *
 *  Analytic source: closed-form experssion
 *
 * ********************************************************************/

#include <unittest++/UnitTest++.h>
#include <lgnSimulator.h>

using namespace lgnSimulator;

void runSystemTest_GRIC(int nt, double dt, int ns, double ds,
                        double C, int wdId, int kxId, int thetaId,
                        const Kernel &W, const Kernel &Kcr,
                        const Kernel &Krg, const Kernel &Kri,const Kernel &Krc,
                        const Kernel &Kig, const Kernel &Kic)
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
    double R0 = 0.0;
    RelayCell relay(integrator, R0);

    //cortical cell
    Interneuron interneuron(integrator);

    //cortical cell
    CorticalCell cortical(integrator);

    //Connect
    relay.addGanglionCell(&ganglion, Krg);
    relay.addCorticalCell(&cortical, Krc);
    relay.addInterNeuron(&interneuron, Kri);
    interneuron.addGanglionCell(&ganglion, Kig);
    interneuron.addCorticalCell(&cortical, Kic);
    cortical.addRelayCell(&relay, Kcr);

    //Compute
    ganglion.computeResponse(&grating);
    relay.computeResponse(&grating);
    interneuron.computeResponse(&grating);
    cortical.computeResponse(&grating);

    cube Rg_e = zeros(r.n_elem, r.n_elem, t.n_elem);
    cube Rr_e = zeros(r.n_elem, r.n_elem, t.n_elem);
    cube Ri_e = zeros(r.n_elem, r.n_elem, t.n_elem);
    cube Rc_e = zeros(r.n_elem, r.n_elem, t.n_elem);

    double kx = spatialFreq*cos(orientation*core::pi/180.);
    double ky = spatialFreq*sin(orientation*core::pi/180.);

    complex<double> Wijl = W.fourierTransform({kx, ky}, wd);
    complex<double> Wr = Wijl*(Krg.fourierTransform({kx, ky}, wd)
                       + Kri.fourierTransform({kx, ky}, wd)
                       * Kig.fourierTransform({kx, ky}, wd))
                       / (1. - Kic.fourierTransform({kx, ky}, wd)
                       * Kcr.fourierTransform({kx, ky}, wd)
                       * Kri.fourierTransform({kx, ky}, wd)
                       - Krc.fourierTransform({kx, ky}, wd)
                       * Kcr.fourierTransform({kx, ky}, wd));
    complex<double> Wc = Wr * Kcr.fourierTransform({kx, ky}, wd);
    complex<double> Wi = Wijl* Kig.fourierTransform({kx, ky}, wd)
            + Wr * Kcr.fourierTransform({kx, ky}, wd)
            * Kic.fourierTransform({kx, ky}, wd);

    for(int l=0; l < int(t.n_elem); l++){
        for(int i = 0; i < int(r.n_elem); i++){
            for(int j = 0; j < int(r.n_elem); j++){
                Rg_e(i,j,l) = C * abs(Wijl)
                        * cos(kx*r[i] + ky*r[j] - wd * t[l] + arg(Wijl));

                Rr_e(i,j,l) = R0 + C * abs(Wr)
                        * cos(kx*r[i] + ky*r[j] - wd * t[l] + arg(Wr));

                Ri_e(i,j,l) = C * abs(Wi)
                        * cos(kx*r[i] + ky*r[j] - wd * t[l] + arg(Wi));

                Rc_e(i,j,l) = C * abs(Wc)
                        * cos(kx*r[i] + ky*r[j] - wd * t[l] + arg(Wc));
            }
        }
    }

    cube diff_g = abs(Rg_e - ganglion.response());
    cube diff_r = abs(Rr_e - relay.response());
    cube diff_i = abs(Ri_e - interneuron.response());
    cube diff_c = abs(Rc_e - cortical.response());

    CHECK_CLOSE(diff_g.max(), 0.0, 1e-10);
    CHECK_CLOSE(diff_r.max(), 0.0, 1e-10);
    CHECK_CLOSE(diff_i.max(), 0.0, 1e-10);
    CHECK_CLOSE(diff_c.max(), 0.0, 1e-10);

}

SUITE(system){

    TEST(runSystemTest_GRIC_0){
        double dt = 0.1;
        double ds = 0.1;
        double w_w = 1.0;
        double w_cr = 1.0;
        double w_rg = 0.5;
        double w_rc = 0.5;
        double w_ri = -0.2;
        double w_ig = 0.5;
        double w_ic = 0.5;

        SpatialGaussian Ws(1, 0.25);
        TemporalDelta Wt(0, dt);
        SeparableKernel W(w_w, &Ws, &Wt);

        SpatialGaussian Kcr_s(1.0, 0.23);
        TemporalDelta Kcr_t(0, dt);
        SeparableKernel Kcr(w_cr, &Kcr_s, &Kcr_t);

        SpatialDelta Krg_s(0.5, ds, {0,0,});
        TemporalDelta Krg_t(0, dt);
        SeparableKernel Krg(w_rg, &Krg_s, &Krg_t);

        SpatialDelta Krc_s(0.5, ds, {0,0,});
        TemporalDelta Krc_t(0, dt);
        SeparableKernel Krc(w_rc, &Krc_s, &Krc_t);

        SpatialGaussian Kri_s(-0.2, 0.23);
        TemporalDelta Kri_t(0, dt);
        SeparableKernel Kri(w_ri, &Kri_s, &Kri_t);

        SpatialDelta Kig_s(0.5, ds, {0,0});
        TemporalDelta Kig_t(0, dt);
        SeparableKernel Kig(w_ig, &Kig_s, &Kig_t);

        SpatialGaussian Kic_s(0.5, 0.23);
        TemporalDelta Kic_t(0, dt);
        SeparableKernel Kic(w_ic, &Kic_s, &Kic_t);

        runSystemTest_GRIC(2, dt, 5, ds,
                           -1, 0, 1, 0,
                           W, Kcr,
                           Krg, Kri, Krc,
                           Kig, Kic);
    }


    TEST(runSystemTest_GRIC_1){
        double dt = 0.1;
        double ds = 0.1;
        double w_w = 1.0;
        double w_cr = -1.5;
        double w_rg = 0.5;
        double w_rc = 0.5;
        double w_ri = -1.5;
        double w_ig = 0.5;
        double w_ic = 1.5;


        SpatialGaussian Ws(1, 0.25);
        TemporalDelta Wt(0, dt);
        SeparableKernel W(w_w, &Ws, &Wt);

        SpatialGaussian Kcr_s(-1.5, 0.23);
        TemporalDelta Kcr_t(0, dt);
        SeparableKernel Kcr(w_cr, &Kcr_s, &Kcr_t);


        SpatialDelta Krg_s(0.5, ds, {0,0,});
        TemporalDelta Krg_t(2*dt, dt);
        SeparableKernel Krg(w_rg, &Krg_s, &Krg_t);

        SpatialDelta Krc_s(0.5, ds, {0,0,});
        TemporalDelta Krc_t(0, dt);
        SeparableKernel Krc(w_rc, &Krc_s, &Krc_t);


        SpatialGaussian Kri_s(-1.5, 0.23);
        TemporalDelta Kri_t(0, dt);
        SeparableKernel Kri(w_ri, &Kri_s, &Kri_t);

        SpatialDelta Kig_s(0.5, ds, {0,0});
        TemporalDelta Kig_t(2*dt, dt);
        SeparableKernel Kig(w_ig, &Kig_s, &Kig_t);

        SpatialGaussian Kic_s(1.5, 0.23);
        TemporalDelta Kic_t(0, dt);
        SeparableKernel Kic(w_ic, &Kic_s, &Kic_t);

        runSystemTest_GRIC(2, dt, 5, ds,
                           -1, 0, 1, 0,
                           W, Kcr,
                           Krg, Kri, Krc,
                           Kig, Kic);
    }

}



