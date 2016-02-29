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

    //DampedOscillator---------------------------------------------------
    TEST(dampedOscillator_test_0) {
        double phaseDuration = 42.5;
        double weight = 0.38;
        double delay = 0.;
        DampedOscillator dampedOsc(phaseDuration, weight, delay);
        CHECK_CLOSE(dampedOsc.temporal(0.5), 0.0369514993891, 1e-10);
        CHECK_CLOSE(real(dampedOsc.fourierTransform(-0.0245436926062)), 22.7988820989, 1e-10);
        CHECK_CLOSE(imag(dampedOsc.fourierTransform(-0.0245436926062)), -3.1173542118, 1e-10);


        CHECK_CLOSE(dampedOsc.temporal(6.0), 0.429120608773, 1e-10);
        CHECK_CLOSE(real(dampedOsc.fourierTransform(-0.294524311274)), -1.1286928585, 1e-10);
        CHECK_CLOSE(imag(dampedOsc.fourierTransform(-0.294524311274)), 0.0062065364, 1e-10);


        CHECK_CLOSE(dampedOsc.temporal(11.5), 0.751331889557, 1e-10);
        CHECK_CLOSE(real(dampedOsc.fourierTransform(-0.564504929942)), -0.3555288674, 1e-10);
        CHECK_CLOSE(imag(dampedOsc.fourierTransform(-0.564504929942)), -0.0651266984, 1e-10);


        CHECK_CLOSE(dampedOsc.temporal(17.0), 0.951056516295, 1e-10);
        CHECK_CLOSE(real(dampedOsc.fourierTransform(-0.83448554861)), -0.0760582312, 1e-10);
        CHECK_CLOSE(imag(dampedOsc.fourierTransform(-0.83448554861)), -0.0917320759, 1e-10);


        CHECK_CLOSE(dampedOsc.temporal(22.5), 0.995734176295, 1e-10);
        CHECK_CLOSE(real(dampedOsc.fourierTransform(-1.10446616728)), -0.0021875183, 1e-10);
        CHECK_CLOSE(imag(dampedOsc.fourierTransform(-1.10446616728)), 0.0152324929, 1e-10);


        CHECK_CLOSE(dampedOsc.temporal(28.0), 0.878081248084, 1e-10);
        CHECK_CLOSE(real(dampedOsc.fourierTransform(-1.37444678595)), -0.0445794299, 1e-10);
        CHECK_CLOSE(imag(dampedOsc.fourierTransform(-1.37444678595)), 0.0315679060, 1e-10);


        CHECK_CLOSE(dampedOsc.temporal(33.5), 0.61727822129, 1e-10);
        CHECK_CLOSE(real(dampedOsc.fourierTransform(-1.64442740461)), -0.0392905913, 1e-10);
        CHECK_CLOSE(imag(dampedOsc.fourierTransform(-1.64442740461)), 0.0014546807, 1e-10);


        CHECK_CLOSE(dampedOsc.temporal(39.0), 0.255842777594, 1e-10);
        CHECK_CLOSE(real(dampedOsc.fourierTransform(-1.91440802328)), -0.0259257859, 1e-10);
        CHECK_CLOSE(imag(dampedOsc.fourierTransform(-1.91440802328)), 0.0006440208, 1e-10);

    }

    TEST(dampedOscillator_test_1) {
        double phaseDuration = 20.6;
        double weight = 0.88;
        double delay = 44.5;
        DampedOscillator dampedOsc(phaseDuration, weight, delay);

        CHECK_CLOSE(dampedOsc.temporal(45.0), 0.0761783767709, 1e-10);
        CHECK_CLOSE(real(dampedOsc.fourierTransform(-2.20893233456)), 0.0355236623, 1e-10);
        CHECK_CLOSE(imag(dampedOsc.fourierTransform(-2.20893233456)), -0.0472345014, 1e-10);


        CHECK_CLOSE(dampedOsc.temporal(50.5), 0.792579042689, 1e-10);
        CHECK_CLOSE(real(dampedOsc.fourierTransform(-2.47891295322)), 0.0327956407, 1e-10);
        CHECK_CLOSE(imag(dampedOsc.fourierTransform(-2.47891295322)), 0.0088992464, 1e-10);


        CHECK_CLOSE(dampedOsc.temporal(56.0), 0.983301195364, 1e-10);
        CHECK_CLOSE(real(dampedOsc.fourierTransform(-2.74889357189)), 0.0044665680, 1e-10);
        CHECK_CLOSE(imag(dampedOsc.fourierTransform(-2.74889357189)), 0.0035356123, 1e-10);


        CHECK_CLOSE(dampedOsc.temporal(61.5), 0.521848255578, 1e-10);
        CHECK_CLOSE(real(dampedOsc.fourierTransform(-3.01887419056)), 0.0192555610, 1e-10);
        CHECK_CLOSE(imag(dampedOsc.fourierTransform(-3.01887419056)), 0.0001958484, 1e-10);


    }

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
