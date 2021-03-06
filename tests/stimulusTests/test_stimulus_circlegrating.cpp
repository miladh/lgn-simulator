/**********************************************************************
 *  Test: Full-field grating
 *
 *  Analytic source: closed-form experssion
 *
 * ********************************************************************/


#include <lgnSimulator.h>
#include <catch.hpp>

using namespace lgnSimulator;



void runStimulusCircleMaskGratingTest(int ns, int nt, double dt, double ds,
                                      double C, int wdId, int kxId, int thetaId,
                                      double maskSize){


    Integrator integrator(nt, dt, ns, ds);
    vec k = integrator.spatialFreqVec();
    vec w = integrator.temporalFreqVec();
    vec orientations = {0., 30., 45., -60., 90., -120., 180., -330.};

    double wd = w(wdId);
    double spatialFreq = k(kxId);
    double orientation = orientations(thetaId);
    double phase =0.0;


    CircleMaskGrating grating(&integrator, spatialFreq, wd, C, phase,
                              orientation, maskSize);
    grating.computeSpatiotemporal();
    grating.computeFourierTransform();

    cx_cube Sdg_ft = integrator.forwardFFT(grating.spatioTemporal());
    cx_cube diff_fourierTransform = grating.fourierTransform() - Sdg_ft;


    cube diff_fourierTransform_real = (real(diff_fourierTransform));
    cube diff_fourierTransform_imag = (imag(diff_fourierTransform));


//    int idc = integrator.nPointsSpatial()/2;
//    cout << "error center: "<<
//            grating.fourierTransform()(idc,idc,0) - Sdg_ft(idc,idc,0)<<endl;

//    cout << diff_fourierTransform_real.max()<< endl;
//    cout << diff_fourierTransform_imag.max()<< endl;

    // Test
    REQUIRE(diff_fourierTransform_real.max() == Approx(0.0).epsilon(1));
    REQUIRE(diff_fourierTransform_imag.max() == Approx(0.0).epsilon(1e-10));

}



TEST_CASE("CircleMaskGrating_test_0"){
     runStimulusCircleMaskGratingTest(9, 1, 1, 0.1,
                                      1.0, 0, 0, 0, 5.0);
}

