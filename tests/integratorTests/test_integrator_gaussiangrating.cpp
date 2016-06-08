/**********************************************************************
 *  Test: 3D inverse fourier transform of full-field grating functions
 *        (cosine waves)
 *
 *  Analytic source: implementation of the spatiotemporal grating
 *                   function in FullFieldGrating class
 *
 * ********************************************************************/

#include <unittest++/UnitTest++.h>
#include <lgnSimulator.h>

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

    GaussianMaskGrating grating(integrator, spatialFreq, orientation,
                                wd, C, maskSize);
    grating.computeFourierTransform();
    grating.computeSpatiotemporal();

    cx_cube grating_fft = integrator.backwardFFT(grating.fourierTransform());
    cx_cube diff = grating.spatioTemporal() - grating_fft;

    cube diff_real = abs(real(diff));
    cube diff_imag = abs(imag(diff));


    uword  idx, idy, idz;
    double maxId = diff_real.max(idx, idy, idz);

    cout << "index: " << idx <<"  "<<  idy << "  "<< idz << endl;

    cout << diff_real.max() << endl;
    cout << diff_imag.max() << endl;

    double idc = integrator.nPointsSpatial()/2;
    cout << "error center: "<<
            grating.spatioTemporal()(idc,idc,0) - real(grating_fft(idc,idc,0))<<endl;

//    cout << grating.spatioTemporal() - real(grating_fft)<<endl;


    //Test
    CHECK_CLOSE(diff_real.max(), 0.0, 1e-9);
    CHECK_CLOSE(diff_imag.max(), 0.0, 1e-9);

}


SUITE(integrator){

    TEST(gaussian_grating_test_0){
        runIntegratorGaussianGratingTest(7, 2, 1, 0.1,
                                         1.0, 0, 0, 0, 0.25);
    }

}
