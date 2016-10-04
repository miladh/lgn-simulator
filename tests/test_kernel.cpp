/**********************************************************************
 *  Test: spatial and temporal kernel functions
 *
 *  Analytic source: by hand
 *
 * ********************************************************************/

#include <lgnSimulator.h>
#include <catch.hpp>

using namespace lgnSimulator;

//Decaying exponential-----------------------------------------



//DifferenceOfExponential--------------------------------------
TEST_CASE("[doe_test_0]") {
    double cenLatency = 0.25;
    double surLatency = 0.125;
    double delay = 0.0;
    DOE doe(cenLatency,surLatency,delay);

    REQUIRE(doe.temporal(-10.24) == Approx( 0.000000000000000).epsilon(1e-14));
    REQUIRE(doe.temporal(0.0) == Approx( 0.000000000000000).epsilon(1e-14));
    REQUIRE(doe.temporal(1.1) == Approx( 0.205469573807283).epsilon(1e-14));
    REQUIRE(doe.temporal(1.8) == Approx( 0.021437459910709).epsilon(1e-14));
    REQUIRE(doe.temporal(2.3) == Approx( 0.003716747228587).epsilon(1e-14));
    REQUIRE(doe.temporal(100.4) == Approx( 0.000000000000000).epsilon(1e-14));

    REQUIRE(real(doe.fourierTransform(-0.0)) == Approx( 0.000000000000000).epsilon(1e-14));
    REQUIRE(imag(doe.fourierTransform(-0.0)) == Approx( 0.000000000000000).epsilon(1e-14));
    REQUIRE(real(doe.fourierTransform(-157.079632679)) == Approx( 0.001926530044585).epsilon(1e-12));
    REQUIRE(imag(doe.fourierTransform(-157.079632679)) == Approx( 0.000229856503266).epsilon(1e-12));
    REQUIRE(real(doe.fourierTransform(314.159265359)) == Approx( 0.000485160381154).epsilon(1e-14));
    REQUIRE(imag(doe.fourierTransform(314.159265359)) == Approx( -0.000028855923396).epsilon(1e-14));
    REQUIRE(real(doe.fourierTransform(0.0383495196971)) == Approx(-0.000206775572248).epsilon(1e-14));
    REQUIRE(imag(doe.fourierTransform(0.0383495196971)) == Approx( 0.009584296015202).epsilon(1e-14));
    REQUIRE(real(doe.fourierTransform(1005.18693069)) == Approx( 0.000047494619442).epsilon(1e-14));
    REQUIRE(imag(doe.fourierTransform(1005.18693069)) == Approx( -0.000000882077205).epsilon(1e-14));

}


TEST_CASE("[doe_test_1]") {
    double cenLatency = 1.43;
    double surLatency = 13.78;
    double delay = 0.4;
    DOE doe(cenLatency,surLatency,delay);


    REQUIRE(doe.temporal(-10.24) == Approx( 0.000000000000000).epsilon(1e-14));
    REQUIRE(doe.temporal(0.0) == Approx( 0.000000000000000).epsilon(1e-14));
    REQUIRE(doe.temporal(1.1) == Approx( 0.206310110952755).epsilon(1e-14));
    REQUIRE(doe.temporal(1.8) == Approx( 0.250540440834636).epsilon(1e-14));
    REQUIRE(doe.temporal(2.3) == Approx( 0.237346187629152).epsilon(1e-14));
    REQUIRE(doe.temporal(100.4) == Approx( -0.000371426190015).epsilon(1e-14));


    REQUIRE(real(doe.fourierTransform(-0.0)) == Approx( 0.000000000000000).epsilon(1e-14));
    REQUIRE(imag(doe.fourierTransform(-0.0)) == Approx( 0.000000000000000).epsilon(1e-14));
    REQUIRE(real(doe.fourierTransform(-157.079632679)) == Approx( -0.000019604682231).epsilon(1e-14));
    REQUIRE(imag(doe.fourierTransform(-157.079632679)) == Approx( -0.000000176262290).epsilon(1e-14));
    REQUIRE(real(doe.fourierTransform(314.159265359)) == Approx( -0.000004901391475).epsilon(1e-14));
    REQUIRE(imag(doe.fourierTransform(314.159265359)) == Approx( 0.000000022033442).epsilon(1e-14));
    REQUIRE(real(doe.fourierTransform(0.0383495196971)) == Approx( 0.558786332919190).epsilon(1e-12));
    REQUIRE(imag(doe.fourierTransform(0.0383495196971)) == Approx( -0.528296608623844).epsilon(1e-12));
    REQUIRE(real(doe.fourierTransform(1005.18693069)) == Approx( -0.000000478165170).epsilon(1e-14));
    REQUIRE(imag(doe.fourierTransform(1005.18693069)) == Approx( 0.000000024164219).epsilon(1e-14));
}


//TwoSidedExponentialDecay--------------------------------------------
TEST_CASE("twoSidedExp_test_0") {
    double tau = 0.5;
    double delay = 9.6;
    TwoSidedExponentialDecay decay(tau,delay);

    REQUIRE(decay.temporal(0.0) ==  Approx(0.000000009174).epsilon(1e-12));
    REQUIRE(real(decay.fourierTransform(-0.0)) ==  Approx(2.0000000000).epsilon(1e-10));
    REQUIRE(imag(decay.fourierTransform(-0.0)) ==  Approx(0.0000000000).epsilon(1e-10));


    REQUIRE(decay.temporal(3.3) ==  Approx(0.000006744030).epsilon(1e-12));
    REQUIRE(real(decay.fourierTransform(-1.79987079112)) ==  Approx(0.0000000000).epsilon(1e-10));
    REQUIRE(imag(decay.fourierTransform(-1.79987079112)) ==  Approx(1.1050433694).epsilon(1e-10));


    REQUIRE(decay.temporal(6.6) ==  Approx(0.004957504353).epsilon(1e-12));
    REQUIRE(real(decay.fourierTransform(-3.59974158224)) ==  Approx(-0.4717498650).epsilon(1e-10));
    REQUIRE(imag(decay.fourierTransform(-3.59974158224)) ==  Approx(0.0000000000).epsilon(1e-10));


    REQUIRE(decay.temporal(9.9) ==  Approx(1.097623272188).epsilon(1e-12));
    REQUIRE(real(decay.fourierTransform(-5.39961237336)) ==  Approx(0.0000000000).epsilon(1e-10));
    REQUIRE(imag(decay.fourierTransform(-5.39961237336)) ==  Approx(-0.2412849841).epsilon(1e-10));


    REQUIRE(decay.temporal(13.2) ==  Approx(0.001493171617).epsilon(1e-12));
    REQUIRE(real(decay.fourierTransform(-7.19948316448)) ==  Approx(0.1432855723).epsilon(1e-10));
    REQUIRE(imag(decay.fourierTransform(-7.19948316448)) ==  Approx(-0.0000000000).epsilon(1e-10));


    REQUIRE(decay.temporal(16.5) ==  Approx(0.000002031263).epsilon(1e-12));
    REQUIRE(real(decay.fourierTransform(-8.9993539556)) ==  Approx(0.0000000000).epsilon(1e-10));
    REQUIRE(imag(decay.fourierTransform(-8.9993539556)) ==  Approx(0.0941305245).epsilon(1e-10));


    REQUIRE(decay.temporal(19.8) ==  Approx(0.000000002763).epsilon(1e-12));
    REQUIRE(real(decay.fourierTransform(10.1447262772)) ==  Approx(-0.0748254664).epsilon(1e-10));
    REQUIRE(imag(decay.fourierTransform(10.1447262772)) ==  Approx(0.0000000000).epsilon(1e-10));


    REQUIRE(decay.temporal(23.1) ==  Approx(0.000000000004).epsilon(1e-12));
    REQUIRE(real(decay.fourierTransform(8.3448554861)) ==  Approx(-0.0000000000).epsilon(1e-10));
    REQUIRE(imag(decay.fourierTransform(8.3448554861)) ==  Approx(-0.1086416073).epsilon(1e-10));

}

TEST_CASE("twoSidedExp_test_1") {
    double tau = 2.24;
    double delay = 7.5;
    TwoSidedExponentialDecay decay(tau,delay);

    REQUIRE(decay.temporal(0.0) ==  Approx(0.015690652100).epsilon(1e-12));
    REQUIRE(real(decay.fourierTransform(-0.0)) ==  Approx(2.0000000000).epsilon(1e-10));
    REQUIRE(imag(decay.fourierTransform(-0.0)) ==  Approx(0.0000000000).epsilon(1e-10));


    REQUIRE(decay.temporal(3.3) ==  Approx(0.068462038770).epsilon(1e-12));
    REQUIRE(real(decay.fourierTransform(-1.79987079112)) ==  Approx(0.0690478125).epsilon(1e-10));
    REQUIRE(imag(decay.fourierTransform(-1.79987079112)) ==  Approx(-0.0931001977).epsilon(1e-10));


    REQUIRE(decay.temporal(6.6) ==  Approx(0.298716122352).epsilon(1e-12));
    REQUIRE(real(decay.fourierTransform(-3.59974158224)) ==  Approx(-0.0087940057).epsilon(1e-10));
    REQUIRE(imag(decay.fourierTransform(-3.59974158224)) ==  Approx(-0.0289899516).epsilon(1e-10));


    REQUIRE(decay.temporal(9.9) ==  Approx(0.152910203167).epsilon(1e-12));
    REQUIRE(real(decay.fourierTransform(-5.39961237336)) ==  Approx(-0.0127847095).epsilon(1e-10));
    REQUIRE(imag(decay.fourierTransform(-5.39961237336)) ==  Approx(-0.0045744422).epsilon(1e-10));


    REQUIRE(decay.temporal(13.2) ==  Approx(0.035045126373).epsilon(1e-12));
    REQUIRE(real(decay.fourierTransform(-7.19948316448)) ==  Approx(-0.0063695833).epsilon(1e-10));
    REQUIRE(imag(decay.fourierTransform(-7.19948316448)) ==  Approx(0.0042560195).epsilon(1e-10));


    REQUIRE(decay.temporal(16.5) ==  Approx(0.008031909298).epsilon(1e-12));
    REQUIRE(real(decay.fourierTransform(-8.9993539556)) ==  Approx(-0.0002409014).epsilon(1e-10));
    REQUIRE(imag(decay.fourierTransform(-8.9993539556)) ==  Approx(0.0049036610).epsilon(1e-10));


    REQUIRE(decay.temporal(19.8) ==  Approx(0.001840814220).epsilon(1e-12));
    REQUIRE(real(decay.fourierTransform(10.1447262772)) ==  Approx(0.0029881229).epsilon(1e-10));
    REQUIRE(imag(decay.fourierTransform(10.1447262772)) ==  Approx(0.0024522891).epsilon(1e-10));


    REQUIRE(decay.temporal(23.1) ==  Approx(0.000421891840).epsilon(1e-12));
    REQUIRE(real(decay.fourierTransform(8.3448554861)) ==  Approx(0.0055365711).epsilon(1e-10));
    REQUIRE(imag(decay.fourierTransform(8.3448554861)) ==  Approx(-0.0013868389).epsilon(1e-10));


    REQUIRE(decay.temporal(26.4) ==  Approx(0.000096692389).epsilon(1e-12));
    REQUIRE(real(decay.fourierTransform(6.54498469498)) ==  Approx(0.0035443816).epsilon(1e-10));
    REQUIRE(imag(decay.fourierTransform(6.54498469498)) ==  Approx(-0.0085568942).epsilon(1e-10));


    REQUIRE(decay.temporal(29.7) ==  Approx(0.000022160699).epsilon(1e-12));
    REQUIRE(real(decay.fourierTransform(4.74511390386)) ==  Approx(-0.0090211812).epsilon(1e-10));
    REQUIRE(imag(decay.fourierTransform(4.74511390386)) ==  Approx(-0.0150509316).epsilon(1e-10));


    REQUIRE(decay.temporal(33.0) ==  Approx(0.000005078958).epsilon(1e-12));
    REQUIRE(real(decay.fourierTransform(2.94524311274)) ==  Approx(-0.0447023346).epsilon(1e-10));
    REQUIRE(imag(decay.fourierTransform(2.94524311274)) ==  Approx(-0.0044027957).epsilon(1e-10));


    REQUIRE(decay.temporal(36.3) ==  Approx(0.000001164034).epsilon(1e-12));
    REQUIRE(real(decay.fourierTransform(1.14537232162)) ==  Approx(-0.1771344329).epsilon(1e-10));
    REQUIRE(imag(decay.fourierTransform(1.14537232162)) ==  Approx(0.1954377296).epsilon(1e-10));


}



//biphasic---------------------------------------------------
TEST_CASE("biphasic_test_0") {
    double phaseDuration = 42.5;
    double weight = 0.38;
    double delay = 0.;
    Biphasic biphasic(phaseDuration,weight,delay);
    REQUIRE(biphasic.temporal(0.5) ==  Approx(0.0369514993891).epsilon(1e-12));
    REQUIRE(real(biphasic.fourierTransform(-0.0245436926062)) ==  Approx(22.7988820989).epsilon(1e-10));
    REQUIRE(imag(biphasic.fourierTransform(-0.0245436926062)) ==  Approx(-3.1173542118).epsilon(1e-10));


    REQUIRE(biphasic.temporal(6.0) ==  Approx(0.429120608773).epsilon(1e-12));
    REQUIRE(real(biphasic.fourierTransform(-0.294524311274)) ==  Approx(-1.1286928585).epsilon(1e-10));
    REQUIRE(imag(biphasic.fourierTransform(-0.294524311274)) ==  Approx(0.0062065364).epsilon(1e-10));


    REQUIRE(biphasic.temporal(11.5) ==  Approx(0.751331889557).epsilon(1e-12));
    REQUIRE(real(biphasic.fourierTransform(-0.564504929942)) ==  Approx(-0.3555288674).epsilon(1e-10));
    REQUIRE(imag(biphasic.fourierTransform(-0.564504929942)) ==  Approx(-0.0651266984).epsilon(1e-10));


    REQUIRE(biphasic.temporal(17.0) ==  Approx(0.951056516295).epsilon(1e-12));
    REQUIRE(real(biphasic.fourierTransform(-0.83448554861)) ==  Approx(-0.0760582312).epsilon(1e-10));
    REQUIRE(imag(biphasic.fourierTransform(-0.83448554861)) ==  Approx(-0.0917320759).epsilon(1e-10));


    REQUIRE(biphasic.temporal(22.5) ==  Approx(0.995734176295).epsilon(1e-12));
    REQUIRE(real(biphasic.fourierTransform(-1.10446616728)) ==  Approx(-0.0021875183).epsilon(1e-10));
    REQUIRE(imag(biphasic.fourierTransform(-1.10446616728)) ==  Approx(0.0152324929).epsilon(1e-10));


    REQUIRE(biphasic.temporal(28.0) ==  Approx(0.878081248084).epsilon(1e-12));
    REQUIRE(real(biphasic.fourierTransform(-1.37444678595)) ==  Approx(-0.0445794299).epsilon(1e-10));
    REQUIRE(imag(biphasic.fourierTransform(-1.37444678595)) ==  Approx(0.0315679060).epsilon(1e-10));


    REQUIRE(biphasic.temporal(33.5) ==  Approx(0.61727822129).epsilon(1e-12));
    REQUIRE(real(biphasic.fourierTransform(-1.64442740461)) ==  Approx(-0.0392905913).epsilon(1e-10));
    REQUIRE(imag(biphasic.fourierTransform(-1.64442740461)) ==  Approx(0.0014546807).epsilon(1e-10));


    REQUIRE(biphasic.temporal(39.0) ==  Approx(0.255842777594).epsilon(1e-12));
    REQUIRE(real(biphasic.fourierTransform(-1.91440802328)) ==  Approx(-0.0259257859).epsilon(1e-10));
    REQUIRE(imag(biphasic.fourierTransform(-1.91440802328)) ==  Approx(0.0006440208).epsilon(1e-10));

}

TEST_CASE("biphasic_test_1") {
    double phaseDuration = 20.6;
    double weight = 0.88;
    double delay = 44.5;
    Biphasic biphasic(phaseDuration,weight,delay);

    REQUIRE(biphasic.temporal(45.0) ==  Approx(0.0761783767709).epsilon(1e-12));
    REQUIRE(real(biphasic.fourierTransform(-2.20893233456)) ==  Approx(0.0355236623).epsilon(1e-10));
    REQUIRE(imag(biphasic.fourierTransform(-2.20893233456)) ==  Approx(-0.0472345014).epsilon(1e-10));


    REQUIRE(biphasic.temporal(50.5) ==  Approx(0.792579042689).epsilon(1e-12));
    REQUIRE(real(biphasic.fourierTransform(-2.47891295322)) ==  Approx(0.0327956407).epsilon(1e-10));
    REQUIRE(imag(biphasic.fourierTransform(-2.47891295322)) ==  Approx(0.0088992464).epsilon(1e-10));


    REQUIRE(biphasic.temporal(56.0) ==  Approx(0.983301195364).epsilon(1e-12));
    REQUIRE(real(biphasic.fourierTransform(-2.74889357189)) ==  Approx(0.0044665680).epsilon(1e-10));
    REQUIRE(imag(biphasic.fourierTransform(-2.74889357189)) ==  Approx(0.0035356123).epsilon(1e-10));


    REQUIRE(biphasic.temporal(61.5) ==  Approx(0.521848255578).epsilon(1e-12));
    REQUIRE(real(biphasic.fourierTransform(-3.01887419056)) ==  Approx(0.0192555610).epsilon(1e-10));
    REQUIRE(imag(biphasic.fourierTransform(-3.01887419056)) ==  Approx(0.0001958484).epsilon(1e-10));


}

//TemporalDelta---------------------------------------------------
TEST_CASE("temporalDelta_test_0") {
    double delay = 1.3;
    TemporalDelta delta(delay,1);
    REQUIRE(delta.temporal(0.5) ==  Approx(0.0).epsilon(1e-12));
    REQUIRE(real(delta.fourierTransform(0.5)) ==  Approx(0.796083798549055).epsilon(1e-12));
    REQUIRE(imag(delta.fourierTransform(0.5)) ==  Approx(0.605186405736039).epsilon(1e-12));

}
TEST_CASE("temporalDelta_test_1") {
    double delay = 0.0;
    TemporalDelta delta(delay,1);
    REQUIRE(delta.temporal(0.0) ==  Approx(1.0).epsilon(1e-12));
    REQUIRE(real(delta.fourierTransform(0.0)) ==  Approx(1.0).epsilon(1e-12));
    REQUIRE(imag(delta.fourierTransform(0.5)) ==  Approx(0.0).epsilon(1e-12));

}
TEST_CASE("temporalDelta_test_2") {
    double delay = -20.3;
    TemporalDelta delta(delay,1);
    REQUIRE(delta.temporal(-20.3) ==  Approx(1.0).epsilon(1e-12));
    REQUIRE(real(delta.fourierTransform(2.3)) ==  Approx(-0.907337331535000).epsilon(1e-12));
    REQUIRE(imag(delta.fourierTransform(2.3)) ==  Approx(-0.420403338239533).epsilon(1e-12));

}

//SpatialDelta---------------------------------------------------
TEST_CASE("spatialDelta_test_0") {
    vec2 shift = {0.0,0.0};
    SpatialDelta delta(1,shift);
    REQUIRE(delta.spatial({0.5, 0.1}) ==  Approx(0.0).epsilon(1e-12));
    REQUIRE(real(delta.fourierTransform({0.5, 0.1})) ==  Approx(1.0).epsilon(1e-12));
    REQUIRE(imag(delta.fourierTransform({0.5, 0.1})) ==  Approx(0.0).epsilon(1e-12));

}
TEST_CASE("spatialDelta_test_1") {
    vec2 shift = {0.5,0.1};
    SpatialDelta delta(1,shift);
    REQUIRE(delta.spatial({0.5, 0.1}) ==  Approx(1).epsilon(1e-12));
    REQUIRE(real(delta.fourierTransform({2.5, -3.1})) == Approx(0.5897880250310923).epsilon(1e-12));
    REQUIRE(imag(delta.fourierTransform({2.5, -3.1})) ==  Approx(-0.8075581004051077).epsilon(1e-12));

}

//Gauss----------------------------------------------------------
TEST_CASE("gaussKernel_test_0") {
    SpatialGaussian G(0.25);
    REQUIRE(G.spatial({0.5,0.1}) ==  Approx(0.079488639761866486).epsilon(1e-12));
    REQUIRE(G.spatial({1.2 ,1.9}) ==  Approx(0).epsilon(1e-12));

    REQUIRE(real(G.fourierTransform({0.5 ,1.1})) ==  Approx(0.9774457376685004).epsilon(1e-12));
    REQUIRE(real(G.fourierTransform({1.5,0.1})) ==  Approx(0.96530371170877705).epsilon(1e-12));

    REQUIRE(imag(G.fourierTransform({0.5,1.1})) == 0.0);
    REQUIRE(imag(G.fourierTransform({1.5,0.1})) == 0.0);
}

//DoG------------------------------------------------------------
TEST_CASE("dogKernel") {
    DOG dog(0.25,0.83,0.85);
    REQUIRE(dog.spatial({0.5,0.1}) ==  Approx(-0.189791527743).epsilon(1e-12));
    REQUIRE(dog.spatial({1.2,1.9}) ==  Approx(-0.00025733892027).epsilon(1e-12));

    REQUIRE(real(dog.fourierTransform({0.5,1.1})) ==  Approx(0.316423256919).epsilon(1e-12));
    REQUIRE(real(dog.fourierTransform({1.5,0.1})) ==  Approx(0.389361200098).epsilon(1e-12));

    REQUIRE(imag(dog.fourierTransform({0.5, 1.1})) == 0.0);
    REQUIRE(imag(dog.fourierTransform({1.5, 0.1})) == 0.0);
}
