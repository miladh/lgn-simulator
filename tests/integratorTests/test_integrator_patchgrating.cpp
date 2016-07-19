/**********************************************************************
 *  Test: 3D inverse fourier transform of patch grating functions
 *
 *  Analytic source: implementation of the spatiotemporal grating
 *                   function in circleMaskGrating class
 *
 * ********************************************************************/


#include <lgnSimulator.h>
#include <catch.hpp>


using namespace lgnSimulator;


void runIntegratorPatchGratingTest(int ns, int nt, double dt, double ds,
                              double C, int wdId, int kxId, int thetaId, double maskSize)
{


    Integrator integrator(nt, dt, ns, ds);
    vec k = integrator.spatialFreqVec();
    vec w = integrator.temporalFreqVec();
    vec orientations = {0., 30., 45., -60., 90., -120., 180., -330.};


    double wd = w(wdId);
    double spatialFreq = k(kxId);
    double orientation = orientations(thetaId);

    CircleMaskGrating grating(&integrator, spatialFreq, orientation, wd, C, maskSize);
    grating.computeFourierTransform();
    grating.computeSpatiotemporal();

    cube grating_fft = integrator.backwardFFT(grating.fourierTransform());
    cube diff = abs(grating.spatioTemporal() - grating_fft);


    int idc = integrator.nPointsSpatial()/2;
    cout << "error center: "<<
            grating.spatioTemporal()(idc,idc,0) - grating_fft(idc,idc,0)<<endl;

    //Test
    REQUIRE(diff.max() ==Approx(0.0).epsilon(1e-4));

}



TEST_CASE("patchGrating_test_0"){
//    runIntegratorPatchGratingTest(9, 1, 1, 0.1,
//                                  1.0, 0, 0, 0, 14.);
}
