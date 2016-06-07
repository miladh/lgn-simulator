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

    //DifferenceOfExponential--------------------------------------
    TEST(doe_test_0) {
        double cenLatency = 0.25;
        double surLatency = 0.125;
        double delay = 0.0;
        DOE doe(cenLatency, surLatency, delay);

        CHECK_CLOSE(doe.temporal(-10.24), 0.000000000000000, 1e-14);
        CHECK_CLOSE(doe.temporal(0.0), 0.000000000000000, 1e-14);
        CHECK_CLOSE(doe.temporal(1.1), 0.205469573807283, 1e-14);
        CHECK_CLOSE(doe.temporal(1.8), 0.021437459910709, 1e-14);
        CHECK_CLOSE(doe.temporal(2.3), 0.003716747228587, 1e-14);
        CHECK_CLOSE(doe.temporal(100.4), 0.000000000000000, 1e-14);

        CHECK_CLOSE(real(doe.fourierTransform(-0.0)), 0.000000000000000, 1e-14);
        CHECK_CLOSE(imag(doe.fourierTransform(-0.0)), 0.000000000000000, 1e-14);
        CHECK_CLOSE(real(doe.fourierTransform(-157.079632679)), 0.001926530044585, 1e-12);
        CHECK_CLOSE(imag(doe.fourierTransform(-157.079632679)), 0.000229856503266, 1e-12);
        CHECK_CLOSE(real(doe.fourierTransform(314.159265359)), 0.000485160381154, 1e-14);
        CHECK_CLOSE(imag(doe.fourierTransform(314.159265359)), -0.000028855923396, 1e-14);
        CHECK_CLOSE(real(doe.fourierTransform(0.0383495196971)),-0.000206775572248, 1e-14);
        CHECK_CLOSE(imag(doe.fourierTransform(0.0383495196971)), 0.009584296015202, 1e-14);
        CHECK_CLOSE(real(doe.fourierTransform(1005.18693069)), 0.000047494619442, 1e-14);
        CHECK_CLOSE(imag(doe.fourierTransform(1005.18693069)), -0.000000882077205, 1e-14);

    }


    TEST(doe_test_1) {
        double cenLatency = 1.43;
        double surLatency = 13.78;
        double delay = 0.4;
        DOE doe(cenLatency, surLatency, delay);


        CHECK_CLOSE(doe.temporal(-10.24), 0.000000000000000, 1e-14);
        CHECK_CLOSE(doe.temporal(0.0), 0.000000000000000, 1e-14);
        CHECK_CLOSE(doe.temporal(1.1), 0.206310110952755, 1e-14);
        CHECK_CLOSE(doe.temporal(1.8), 0.250540440834636, 1e-14);
        CHECK_CLOSE(doe.temporal(2.3), 0.237346187629152, 1e-14);
        CHECK_CLOSE(doe.temporal(100.4), -0.000371426190015, 1e-14);


        CHECK_CLOSE(real(doe.fourierTransform(-0.0)), 0.000000000000000, 1e-14);
        CHECK_CLOSE(imag(doe.fourierTransform(-0.0)), 0.000000000000000, 1e-14);
        CHECK_CLOSE(real(doe.fourierTransform(-157.079632679)), -0.000019604682231, 1e-14);
        CHECK_CLOSE(imag(doe.fourierTransform(-157.079632679)), -0.000000176262290, 1e-14);
        CHECK_CLOSE(real(doe.fourierTransform(314.159265359)), -0.000004901391475, 1e-14);
        CHECK_CLOSE(imag(doe.fourierTransform(314.159265359)), 0.000000022033442, 1e-14);
        CHECK_CLOSE(real(doe.fourierTransform(0.0383495196971)), 0.558786332919190, 1e-12);
        CHECK_CLOSE(imag(doe.fourierTransform(0.0383495196971)),-0.528296608623844, 1e-12);
        CHECK_CLOSE(real(doe.fourierTransform(1005.18693069)), -0.000000478165170, 1e-14);
        CHECK_CLOSE(imag(doe.fourierTransform(1005.18693069)), 0.000000024164219, 1e-14);


    }


    //TwoSidedExponentialDecay--------------------------------------------
    TEST(twoSidedExp_test_0) {
        double tau = 0.5;
        double delay = 9.6;
        TwoSidedExponentialDecay decay(tau, delay);

        CHECK_CLOSE(decay.temporal(0.0), 0.000000009174, 1e-12);
        CHECK_CLOSE(real(decay.fourierTransform(-0.0)), 2.0000000000, 1e-10);
        CHECK_CLOSE(imag(decay.fourierTransform(-0.0)), 0.0000000000, 1e-10);


        CHECK_CLOSE(decay.temporal(3.3), 0.000006744030, 1e-12);
        CHECK_CLOSE(real(decay.fourierTransform(-1.79987079112)), 0.0000000000, 1e-10);
        CHECK_CLOSE(imag(decay.fourierTransform(-1.79987079112)), 1.1050433694, 1e-10);


        CHECK_CLOSE(decay.temporal(6.6), 0.004957504353, 1e-12);
        CHECK_CLOSE(real(decay.fourierTransform(-3.59974158224)), -0.4717498650, 1e-10);
        CHECK_CLOSE(imag(decay.fourierTransform(-3.59974158224)), 0.0000000000, 1e-10);


        CHECK_CLOSE(decay.temporal(9.9), 1.097623272188, 1e-12);
        CHECK_CLOSE(real(decay.fourierTransform(-5.39961237336)), 0.0000000000, 1e-10);
        CHECK_CLOSE(imag(decay.fourierTransform(-5.39961237336)), -0.2412849841, 1e-10);


        CHECK_CLOSE(decay.temporal(13.2), 0.001493171617, 1e-12);
        CHECK_CLOSE(real(decay.fourierTransform(-7.19948316448)), 0.1432855723, 1e-10);
        CHECK_CLOSE(imag(decay.fourierTransform(-7.19948316448)), -0.0000000000, 1e-10);


        CHECK_CLOSE(decay.temporal(16.5), 0.000002031263, 1e-12);
        CHECK_CLOSE(real(decay.fourierTransform(-8.9993539556)), 0.0000000000, 1e-10);
        CHECK_CLOSE(imag(decay.fourierTransform(-8.9993539556)), 0.0941305245, 1e-10);


        CHECK_CLOSE(decay.temporal(19.8), 0.000000002763, 1e-12);
        CHECK_CLOSE(real(decay.fourierTransform(10.1447262772)), -0.0748254664, 1e-10);
        CHECK_CLOSE(imag(decay.fourierTransform(10.1447262772)), 0.0000000000, 1e-10);


        CHECK_CLOSE(decay.temporal(23.1), 0.000000000004, 1e-12);
        CHECK_CLOSE(real(decay.fourierTransform(8.3448554861)), -0.0000000000, 1e-10);
        CHECK_CLOSE(imag(decay.fourierTransform(8.3448554861)), -0.1086416073, 1e-10);

    }

    TEST(twoSidedExp_test_1) {
        double tau = 2.24;
        double delay = 7.5;
        TwoSidedExponentialDecay decay(tau, delay);

        CHECK_CLOSE(decay.temporal(0.0), 0.015690652100, 1e-12);
        CHECK_CLOSE(real(decay.fourierTransform(-0.0)), 2.0000000000, 1e-10);
        CHECK_CLOSE(imag(decay.fourierTransform(-0.0)), 0.0000000000, 1e-10);


        CHECK_CLOSE(decay.temporal(3.3), 0.068462038770, 1e-12);
        CHECK_CLOSE(real(decay.fourierTransform(-1.79987079112)), 0.0690478125, 1e-10);
        CHECK_CLOSE(imag(decay.fourierTransform(-1.79987079112)), -0.0931001977, 1e-10);


        CHECK_CLOSE(decay.temporal(6.6), 0.298716122352, 1e-12);
        CHECK_CLOSE(real(decay.fourierTransform(-3.59974158224)), -0.0087940057, 1e-10);
        CHECK_CLOSE(imag(decay.fourierTransform(-3.59974158224)), -0.0289899516, 1e-10);


        CHECK_CLOSE(decay.temporal(9.9), 0.152910203167, 1e-12);
        CHECK_CLOSE(real(decay.fourierTransform(-5.39961237336)), -0.0127847095, 1e-10);
        CHECK_CLOSE(imag(decay.fourierTransform(-5.39961237336)), -0.0045744422, 1e-10);


        CHECK_CLOSE(decay.temporal(13.2), 0.035045126373, 1e-12);
        CHECK_CLOSE(real(decay.fourierTransform(-7.19948316448)), -0.0063695833, 1e-10);
        CHECK_CLOSE(imag(decay.fourierTransform(-7.19948316448)), 0.0042560195, 1e-10);


        CHECK_CLOSE(decay.temporal(16.5), 0.008031909298, 1e-12);
        CHECK_CLOSE(real(decay.fourierTransform(-8.9993539556)), -0.0002409014, 1e-10);
        CHECK_CLOSE(imag(decay.fourierTransform(-8.9993539556)), 0.0049036610, 1e-10);


        CHECK_CLOSE(decay.temporal(19.8), 0.001840814220, 1e-12);
        CHECK_CLOSE(real(decay.fourierTransform(10.1447262772)), 0.0029881229, 1e-10);
        CHECK_CLOSE(imag(decay.fourierTransform(10.1447262772)), 0.0024522891, 1e-10);


        CHECK_CLOSE(decay.temporal(23.1), 0.000421891840, 1e-12);
        CHECK_CLOSE(real(decay.fourierTransform(8.3448554861)), 0.0055365711, 1e-10);
        CHECK_CLOSE(imag(decay.fourierTransform(8.3448554861)), -0.0013868389, 1e-10);


        CHECK_CLOSE(decay.temporal(26.4), 0.000096692389, 1e-12);
        CHECK_CLOSE(real(decay.fourierTransform(6.54498469498)), 0.0035443816, 1e-10);
        CHECK_CLOSE(imag(decay.fourierTransform(6.54498469498)), -0.0085568942, 1e-10);


        CHECK_CLOSE(decay.temporal(29.7), 0.000022160699, 1e-12);
        CHECK_CLOSE(real(decay.fourierTransform(4.74511390386)), -0.0090211812, 1e-10);
        CHECK_CLOSE(imag(decay.fourierTransform(4.74511390386)), -0.0150509316, 1e-10);


        CHECK_CLOSE(decay.temporal(33.0), 0.000005078958, 1e-12);
        CHECK_CLOSE(real(decay.fourierTransform(2.94524311274)), -0.0447023346, 1e-10);
        CHECK_CLOSE(imag(decay.fourierTransform(2.94524311274)), -0.0044027957, 1e-10);


        CHECK_CLOSE(decay.temporal(36.3), 0.000001164034, 1e-12);
        CHECK_CLOSE(real(decay.fourierTransform(1.14537232162)), -0.1771344329, 1e-10);
        CHECK_CLOSE(imag(decay.fourierTransform(1.14537232162)), 0.1954377296, 1e-10);


    }



    //biphasic---------------------------------------------------
    TEST(biphasic_test_0) {
        double phaseDuration = 42.5;
        double weight = 0.38;
        double delay = 0.;
        Biphasic biphasic(phaseDuration, weight, delay);
        CHECK_CLOSE(biphasic.temporal(0.5), 0.0369514993891, 1e-12);
        CHECK_CLOSE(real(biphasic.fourierTransform(-0.0245436926062)), 22.7988820989, 1e-10);
        CHECK_CLOSE(imag(biphasic.fourierTransform(-0.0245436926062)), -3.1173542118, 1e-10);


        CHECK_CLOSE(biphasic.temporal(6.0), 0.429120608773, 1e-12);
        CHECK_CLOSE(real(biphasic.fourierTransform(-0.294524311274)), -1.1286928585, 1e-10);
        CHECK_CLOSE(imag(biphasic.fourierTransform(-0.294524311274)), 0.0062065364, 1e-10);


        CHECK_CLOSE(biphasic.temporal(11.5), 0.751331889557, 1e-12);
        CHECK_CLOSE(real(biphasic.fourierTransform(-0.564504929942)), -0.3555288674, 1e-10);
        CHECK_CLOSE(imag(biphasic.fourierTransform(-0.564504929942)), -0.0651266984, 1e-10);


        CHECK_CLOSE(biphasic.temporal(17.0), 0.951056516295, 1e-12);
        CHECK_CLOSE(real(biphasic.fourierTransform(-0.83448554861)), -0.0760582312, 1e-10);
        CHECK_CLOSE(imag(biphasic.fourierTransform(-0.83448554861)), -0.0917320759, 1e-10);


        CHECK_CLOSE(biphasic.temporal(22.5), 0.995734176295, 1e-12);
        CHECK_CLOSE(real(biphasic.fourierTransform(-1.10446616728)), -0.0021875183, 1e-10);
        CHECK_CLOSE(imag(biphasic.fourierTransform(-1.10446616728)), 0.0152324929, 1e-10);


        CHECK_CLOSE(biphasic.temporal(28.0), 0.878081248084, 1e-12);
        CHECK_CLOSE(real(biphasic.fourierTransform(-1.37444678595)), -0.0445794299, 1e-10);
        CHECK_CLOSE(imag(biphasic.fourierTransform(-1.37444678595)), 0.0315679060, 1e-10);


        CHECK_CLOSE(biphasic.temporal(33.5), 0.61727822129, 1e-12);
        CHECK_CLOSE(real(biphasic.fourierTransform(-1.64442740461)), -0.0392905913, 1e-10);
        CHECK_CLOSE(imag(biphasic.fourierTransform(-1.64442740461)), 0.0014546807, 1e-10);


        CHECK_CLOSE(biphasic.temporal(39.0), 0.255842777594, 1e-12);
        CHECK_CLOSE(real(biphasic.fourierTransform(-1.91440802328)), -0.0259257859, 1e-10);
        CHECK_CLOSE(imag(biphasic.fourierTransform(-1.91440802328)), 0.0006440208, 1e-10);

    }

    TEST(biphasic_test_1) {
        double phaseDuration = 20.6;
        double weight = 0.88;
        double delay = 44.5;
        Biphasic biphasic(phaseDuration, weight, delay);

        CHECK_CLOSE(biphasic.temporal(45.0), 0.0761783767709, 1e-12);
        CHECK_CLOSE(real(biphasic.fourierTransform(-2.20893233456)), 0.0355236623, 1e-10);
        CHECK_CLOSE(imag(biphasic.fourierTransform(-2.20893233456)), -0.0472345014, 1e-10);


        CHECK_CLOSE(biphasic.temporal(50.5), 0.792579042689, 1e-12);
        CHECK_CLOSE(real(biphasic.fourierTransform(-2.47891295322)), 0.0327956407, 1e-10);
        CHECK_CLOSE(imag(biphasic.fourierTransform(-2.47891295322)), 0.0088992464, 1e-10);


        CHECK_CLOSE(biphasic.temporal(56.0), 0.983301195364, 1e-12);
        CHECK_CLOSE(real(biphasic.fourierTransform(-2.74889357189)), 0.0044665680, 1e-10);
        CHECK_CLOSE(imag(biphasic.fourierTransform(-2.74889357189)), 0.0035356123, 1e-10);


        CHECK_CLOSE(biphasic.temporal(61.5), 0.521848255578, 1e-12);
        CHECK_CLOSE(real(biphasic.fourierTransform(-3.01887419056)), 0.0192555610, 1e-10);
        CHECK_CLOSE(imag(biphasic.fourierTransform(-3.01887419056)), 0.0001958484, 1e-10);


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
        vec2 shift = {0.0, 0.0};
        SpatialDelta delta(1, shift);
        CHECK_CLOSE(delta.spatial({0.5, 0.1}), 0.0, 1e-12);
        CHECK_CLOSE(real(delta.fourierTransform({0.5, 0.1})), 1.0, 1e-12);
        CHECK_CLOSE(imag(delta.fourierTransform({0.5, 0.1})), 0.0, 1e-12);

    }
    TEST(spatialDelta_test_1) {
        //w = -1.3
        vec2 shift = {0.5, 0.1};
        SpatialDelta delta(1, shift);
        CHECK_CLOSE(delta.spatial({0.5, 0.1}), 1, 1e-12);
        CHECK_CLOSE(real(delta.fourierTransform({2.5, -3.1})),0.5897880250310923, 1e-12);
        CHECK_CLOSE(imag(delta.fourierTransform({2.5, -3.1})), -0.8075581004051077,1e-12);

    }

    //Gauss----------------------------------------------------------
    TEST(gaussKernel_test_0) {
        SpatialGaussian G(0.25);
        CHECK_CLOSE(G.spatial({0.5, 0.1}), 0.079488639761866486, 1e-12);
        CHECK_CLOSE(G.spatial({1.2, 1.9}), 0, 1e-12);

        CHECK_CLOSE(real(G.fourierTransform({0.5, 1.1})), 0.9774457376685004, 1e-12);
        CHECK_CLOSE(real(G.fourierTransform({1.5, 0.1})), 0.96530371170877705, 1e-12);

        CHECK_EQUAL(imag(G.fourierTransform({0.5, 1.1})), 0.0);
        CHECK_EQUAL(imag(G.fourierTransform({1.5, 0.1})), 0.0);
    }
    TEST(gaussKernel_test_1) {
        SpatialGaussian G(0.25);
        CHECK_CLOSE(G.spatial({0.5, 0.1}),0.07948863976186649, 1e-12);
        CHECK_CLOSE(G.spatial({1.2, 1.9}), 0, 1e-12);

        CHECK_CLOSE(real(G.fourierTransform({0.5, 1.1})),0.9774457376685004, 1e-12);
        CHECK_CLOSE(real(G.fourierTransform({1.5, 0.1})),0.9653037117087772, 1e-12);

        CHECK_EQUAL(imag(G.fourierTransform({0.5, 1.1})), 0.0);
        CHECK_EQUAL(imag(G.fourierTransform({1.5, 0.1})), 0.0);
    }

    //DoG------------------------------------------------------------
    TEST(dogKernel) {
        DOG dog(0.25, 0.83, 0.85);
        CHECK_CLOSE(dog.spatial({0.5, 0.1}), -0.189791527743, 1e-12);
        CHECK_CLOSE(dog.spatial({1.2, 1.9}), -0.00025733892027, 1e-12);

        CHECK_CLOSE(real(dog.fourierTransform({0.5, 1.1})), 0.316423256919, 1e-12);
        CHECK_CLOSE(real(dog.fourierTransform({1.5, 0.1})), 0.389361200098, 1e-12);

        CHECK_EQUAL(imag(dog.fourierTransform({0.5, 1.1})), 0.0);
        CHECK_EQUAL(imag(dog.fourierTransform({1.5, 0.1})), 0.0);
    }

}
