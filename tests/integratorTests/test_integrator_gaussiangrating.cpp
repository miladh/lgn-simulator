/**********************************************************************
 *  Test: 3D inverse fourier transform of gaussian patch grating f
 *        unctions
 *
 *  Analytic source: implementation of the spatiotemporal grating
 *                   function in GaussianMaskGrating class
 *
 * ********************************************************************/
#include <lgnSimulator.h>
#include <catch.hpp>

using namespace lgnSimulator;


void runIntegratorGaussianGratingTest(int ns, int nt, double dt, double ds,
                                      double C, int wdId, int kxId, int thetaId,
                                      double maskSize)
{


    Integrator integrator(nt, dt, ns, ds);
    vec k = integrator.spatialFreqVec();
    vec w = integrator.temporalFreqVec();
    vec orientations = {0., 30., 45., -60., 90., -120., 180., -330.};


    double wd = w(wdId);
    double spatialFreq = k(kxId);
    double orientation = orientations(thetaId);
    double phase = 0.0;

    GaussianMaskGrating grating(&integrator, spatialFreq, wd, C, phase,
                                orientation, maskSize);
    grating.computeFourierTransform();
    grating.computeSpatiotemporal();

    cube grating_fft = integrator.backwardFFT(grating.fourierTransform());
    cube diff = grating.spatioTemporal() - grating_fft;

    cube diff_real = abs(real(diff));
    cube diff_imag = abs(imag(diff));


    //    cout << diff_real.max() << endl;
    //    cout << diff_imag.max() << endl;

    //Test
    REQUIRE(diff_real.max() == Approx(0.0).epsilon(1e-5));
    REQUIRE(diff_imag.max() == Approx(0.0).epsilon(1e-5));

}



TEST_CASE("gaussian_grating_test_0"){
    runIntegratorGaussianGratingTest(7, 2, 1, 0.05,
                                     1.0, 0, 0, 0, 0.2);
}


TEST_CASE("gaussian_grating_test_1"){
    runIntegratorGaussianGratingTest(6, 2, 1, 0.1,
                                     1.0, 0, 0, 0, 0.3);
}


TEST_CASE("gaussian_grating_test_2"){
    runIntegratorGaussianGratingTest(6, 2, 1, 0.1,
                                     1.0, 0, 0, 0, 0.4);
}


TEST_CASE("gaussian_grating_test_3"){
    runIntegratorGaussianGratingTest(6, 2, 1, 0.1,
                                     1.0, 0, 0, 0, 0.5);
}


TEST_CASE("gaussian_grating_test_4"){
    runIntegratorGaussianGratingTest(6, 2, 1, 0.1,
                                     1.0, 0, 0, 0, 0.7);
}


TEST_CASE("gaussian_grating_test_5"){
    runIntegratorGaussianGratingTest(6, 2, 1, 0.1,
                                     1.0, 0, 0, 0, 0.9);
}



TEST_CASE("gaussian_grating_test_6"){
    runIntegratorGaussianGratingTest(7, 2, 1, 0.05,
                                     1.0, 1, 0, 0, 0.2);
}


TEST_CASE("gaussian_grating_test_7"){
    runIntegratorGaussianGratingTest(6, 3, 1, 0.1,
                                     1.0, 2, 0, 0, 0.3);
}


TEST_CASE("gaussian_grating_test_8"){
    runIntegratorGaussianGratingTest(6, 2, 1, 0.1,
                                     1.0, 3, 0, 0, 0.4);
}


TEST_CASE("gaussian_grating_test_9"){
    runIntegratorGaussianGratingTest(6, 2, 1, 0.1,
                                     1.0, 3, 0, 0, 0.5);
}












































