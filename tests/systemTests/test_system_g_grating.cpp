/**********************************************************************
 *  Test: response of ganglion cell with full-field grating stimulus
 *        R(r,t) = W(r,t) * S(r,t)
 *
 *  Analytic source: closed-form experssion
 *
 * ********************************************************************/

#include <lgnSimulator.h>
#include <catch.hpp>


using namespace lgnSimulator;


void runSystemTest_G(int nt, double dt, int ns, double ds,
                     double C, int wdId, int kxId, int thetaId,
                     const Kernel &W)
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
    double phase = orientation+ 18.;

    FullFieldGrating grating(&integrator, spatialFreq, wd, C, phase, orientation);
    grating.computeFourierTransform();

    //ganglion cell
    GanglionCell ganglion(&integrator, W);
    ganglion.computeResponse(&grating);

    cube Rg_e = zeros(r.n_elem, r.n_elem, t.n_elem);

    double kx = spatialFreq*cos(orientation*core::pi/180.);
    double ky = spatialFreq*sin(orientation*core::pi/180.);
    complex<double> Wijl = W.fourierTransform({kx, ky}, wd);

    for(int l=0; l < int(t.n_elem); l++){
        for(int i = 0; i < int(r.n_elem); i++){
            for(int j = 0; j < int(r.n_elem); j++){
                Rg_e(i,j,l) = C * abs(Wijl)
                        * cos(kx*r[i] + ky*r[j] - wd * t[l] + arg(Wijl)
                              + phase*core::pi/180.);
            }
        }
    }

    cube diff = abs(Rg_e - ganglion.response());
    REQUIRE(diff.max()== Approx(0.0).epsilon(1e-10));

}



TEST_CASE("runSystemTest_G_0"){
    double dt = 0.1;
    double weight = 1.0;
    SpatialGaussian Ws(0.25);
    TemporalDelta Wt(0, dt);
    SeparableKernel W(weight, &Ws, &Wt);

    runSystemTest_G(2, dt, 2, 0.1,
                    -1, 0, 1, 0, W);
}


TEST_CASE("runSystemTest_G_1"){
    double dt = 0.1;
    double weight = 1.0;
    SpatialGaussian Ws( 0.25);
    TemporalDelta Wt(1*dt, dt);
    SeparableKernel W(weight, &Ws, &Wt);

    runSystemTest_G(3, dt, 2, 0.1,
                    -1, 3, 1, 4, W);
}







