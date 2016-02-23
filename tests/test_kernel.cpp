/**********************************************************************
 *  Test: spatial and temporal kernel functions
 *
 *  Analytic source: by hand
 *
 * ********************************************************************/

#include <unittest++/UnitTest++.h>
#include <lgnSimulator.h>

using namespace lgnSimulator;

SUITE(kernel){

    //TemporalDelta---------------------------------------------------
    TEST(temporalDelta_test_0) {
        double delay = 1.3;
        TemporalDelta delta(delay, 1);
        CHECK_CLOSE(delta.temporal(0.5), 0.0, 1e-12);
        CHECK_CLOSE(real(delta.fourierTransform(0.5)),0.796083798549055, 1e-12);
        CHECK_CLOSE(imag(delta.fourierTransform(0.5)),0.605186405736039, 1e-12);

    }
    TEST(temporalDelta_test_1) {
        double delay = 0.0;
        TemporalDelta delta(delay, 1);
        CHECK_CLOSE(delta.temporal(0.0), 1.0, 1e-12);
        CHECK_CLOSE(real(delta.fourierTransform(0.0)),1.0, 1e-12);
        CHECK_CLOSE(imag(delta.fourierTransform(0.5)),0.0, 1e-12);

    }
    TEST(temporalDelta_test_2) {
        double delay = -20.3;
        TemporalDelta delta(delay, 1);
        CHECK_CLOSE(delta.temporal(-20.3), 1.0, 1e-12);
        CHECK_CLOSE(real(delta.fourierTransform(2.3)),-0.907337331535000, 1e-12);
        CHECK_CLOSE(imag(delta.fourierTransform(2.3)),-0.420403338239533, 1e-12);

    }

    //SpatialDelta---------------------------------------------------
    TEST(spatialDelta_test_0) {
        double w = 1.3;
        vec2 shift = {0.0, 0.0};
        SpatialDelta delta(w, 1, shift);
        CHECK_CLOSE(delta.spatial({0.5, 0.1}), 0.0, 1e-12);
        CHECK_CLOSE(real(delta.fourierTransform({0.5, 0.1})), w, 1e-12);
        CHECK_CLOSE(imag(delta.fourierTransform({0.5, 0.1})), 0.0, 1e-12);

    }
    TEST(spatialDelta_test_1) {
        double w = -1.3;
        vec2 shift = {0.5, 0.1};
        SpatialDelta delta(w,1, shift);
        CHECK_CLOSE(delta.spatial({0.5, 0.1}), w, 1e-12);
        CHECK_CLOSE(real(delta.fourierTransform({2.5, -3.1})),-0.76672443254042, 1e-12);
        CHECK_CLOSE(imag(delta.fourierTransform({2.5, -3.1})), 1.04982553052664, 1e-12);

    }

    //Gauss----------------------------------------------------------
    TEST(gaussKernel_test_0) {
        Gaussian G(1.0, 0.25);
        CHECK_CLOSE(G.spatial({0.5, 0.1}), 0.079488639761866486, 1e-12);
        CHECK_CLOSE(G.spatial({1.2, 1.9}), 0, 1e-12);

        CHECK_CLOSE(real(G.fourierTransform({0.5, 1.1})), 0.9774457376685004, 1e-12);
        CHECK_CLOSE(real(G.fourierTransform({1.5, 0.1})), 0.96530371170877705, 1e-12);

        CHECK_EQUAL(imag(G.fourierTransform({0.5, 1.1})), 0.0);
        CHECK_EQUAL(imag(G.fourierTransform({1.5, 0.1})), 0.0);
    }
    TEST(gaussKernel_test_1) {
        Gaussian G(-0.75, 0.25);
        CHECK_CLOSE(G.spatial({0.5, 0.1}), -0.059616479821399865, 1e-12);
        CHECK_CLOSE(G.spatial({1.2, 1.9}), 0, 1e-12);

        CHECK_CLOSE(real(G.fourierTransform({0.5, 1.1})), -0.7330843032513753, 1e-12);
        CHECK_CLOSE(real(G.fourierTransform({1.5, 0.1})), -0.72397778378158284, 1e-12);

        CHECK_EQUAL(imag(G.fourierTransform({0.5, 1.1})), 0.0);
        CHECK_EQUAL(imag(G.fourierTransform({1.5, 0.1})), 0.0);
    }

    //DoG------------------------------------------------------------
    TEST(dogKernel) {
        DOG dog(1.0, 0.25, 0.85, 0.83);
        CHECK_CLOSE(dog.spatial({0.5, 0.1}), -0.189791527743, 1e-12);
        CHECK_CLOSE(dog.spatial({1.2, 1.9}), -0.00025733892027, 1e-12);

        CHECK_CLOSE(real(dog.fourierTransform({0.5, 1.1})), 0.316423256919, 1e-12);
        CHECK_CLOSE(real(dog.fourierTransform({1.5, 0.1})), 0.389361200098, 1e-12);

        CHECK_EQUAL(imag(dog.fourierTransform({0.5, 1.1})), 0.0);
        CHECK_EQUAL(imag(dog.fourierTransform({1.5, 0.1})), 0.0);
    }

}
