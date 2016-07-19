/**********************************************************************
 *  Test: 3D inverse fourier transform of a constant functions
 *        F(r,t) = F(r) * F(t) = C * C
 *
 *  Analytic source: closed-form experssion
 *
 * ********************************************************************/

#include <lgnSimulator.h>
#include <catch.hpp>

using namespace lgnSimulator;

void runTest(int ns, int nt, double dt, double ds)
{

    Integrator integrator(nt, dt, ns, ds);
    vec r = integrator.spatialVec();
    vec k = integrator.spatialFreqVec();
    vec t = integrator.timeVec();
    vec w = integrator.temporalFreqVec();

    TemporallyConstant Ft(integrator.timeInterval(),
                          integrator.temporalFreqResolution());
    SpatiallyConstant  Fs(integrator.lengthInterval(), integrator.spatialFreqResolution());

    cube F = zeros<cube>(k.n_elem, k.n_elem, w.n_elem);
    cx_cube G = zeros<cx_cube>(k.n_elem, k.n_elem, w.n_elem);

    for(int l = 0; l < int(w.n_elem); l++){
        for(int i = 0; i < int(k.n_elem); i++){
            for(int j = 0; j < int(k.n_elem); j++){
                F(i,j,l) = Fs.spatial(vec2{r[i], r[j]})
                        * Ft.temporal(t[l]);
                G(i,j,l) = Fs.fourierTransform({k[i], k[j]})
                        * Ft.fourierTransform(w[l]);
            }
        }
    }

    // Backward
    cube F_fft = integrator.backwardFFT(G);
    cube diff = F - F_fft;

    cube diff_real = abs(real(diff));
    cube diff_imag = abs(imag(diff));


    //Test
    REQUIRE(diff_real.max() == Approx(0.0).epsilon(1e-9));
    REQUIRE(diff_imag.max() == Approx(0.0).epsilon(1e-9));


}


TEST_CASE("constant_0"){
    runTest(2, 1, 0.01, 0.355);
}

TEST_CASE("constant_1"){
    runTest(2, 1, 0.255, 0.015);
}

TEST_CASE("constant_2"){
    runTest(2, 1, 0.5, 0.2);
}

TEST_CASE("constant_3"){
    runTest(2, 2, 0.01, 0.355);
}

TEST_CASE("constant_4"){
    runTest(2, 2, 0.255, 0.015);
}

TEST_CASE("constant_5"){
    runTest(2, 2, 0.5, 0.2);
}

TEST_CASE("constant_6"){
    runTest(3, 1, 0.01, 0.355);
}

TEST_CASE("constant_7"){
    runTest(3, 1, 0.255, 0.015);
}

TEST_CASE("constant_8"){
    runTest(3, 1, 0.5, 0.2);
}

TEST_CASE("constant_9"){
    runTest(3, 2, 0.01, 0.355);
}

TEST_CASE("constant_10"){
    runTest(3, 2, 0.255, 0.015);
}

TEST_CASE("constant_11"){
    runTest(3, 2, 0.5, 0.2);
}

TEST_CASE("constant_12"){
    runTest(2, 1, 0.01, 0.355);
}

TEST_CASE("constant_13"){
    runTest(2, 1, 0.255, 0.015);
}

TEST_CASE("constant_14"){
    runTest(2, 1, 0.5, 0.2);
}

TEST_CASE("constant_15"){
    runTest(2, 2, 0.01, 0.355);
}

TEST_CASE("constant_16"){
    runTest(2, 2, 0.255, 0.015);
}

TEST_CASE("constant_17"){
    runTest(2, 2, 0.5, 0.2);
}

TEST_CASE("constant_18"){
    runTest(3, 1, 0.01, 0.355);
}

TEST_CASE("constant_19"){
    runTest(3, 1, 0.255, 0.015);
}

TEST_CASE("constant_20"){
    runTest(3, 1, 0.5, 0.2);
}

TEST_CASE("constant_21"){
    runTest(3, 2, 0.01, 0.355);
}

TEST_CASE("constant_22"){
    runTest(3, 2, 0.255, 0.015);
}

TEST_CASE("constant_23"){
    runTest(3, 2, 0.5, 0.2);
}

TEST_CASE("constant_24"){
    runTest(2, 1, 0.01, 0.355);
}

TEST_CASE("constant_25"){
    runTest(2, 1, 0.255, 0.015);
}

TEST_CASE("constant_26"){
    runTest(2, 1, 0.5, 0.2);
}

TEST_CASE("constant_27"){
    runTest(2, 2, 0.01, 0.355);
}

TEST_CASE("constant_28"){
    runTest(2, 2, 0.255, 0.015);
}

TEST_CASE("constant_29"){
    runTest(2, 2, 0.5, 0.2);
}

TEST_CASE("constant_30"){
    runTest(3, 1, 0.01, 0.355);
}

TEST_CASE("constant_31"){
    runTest(3, 1, 0.255, 0.015);
}

TEST_CASE("constant_32"){
    runTest(3, 1, 0.5, 0.2);
}

TEST_CASE("constant_33"){
    runTest(3, 2, 0.01, 0.355);
}

TEST_CASE("constant_34"){
    runTest(3, 2, 0.255, 0.015);
}

TEST_CASE("constant_35"){
    runTest(3, 2, 0.5, 0.2);
}

TEST_CASE("constant_36"){
    runTest(2, 1, 0.01, 0.355);
}

TEST_CASE("constant_37"){
    runTest(2, 1, 0.255, 0.015);
}

TEST_CASE("constant_38"){
    runTest(2, 1, 0.5, 0.2);
}

TEST_CASE("constant_39"){
    runTest(2, 2, 0.01, 0.355);
}

TEST_CASE("constant_40"){
    runTest(2, 2, 0.255, 0.015);
}

TEST_CASE("constant_41"){
    runTest(2, 2, 0.5, 0.2);
}

TEST_CASE("constant_42"){
    runTest(3, 1, 0.01, 0.355);
}

TEST_CASE("constant_43"){
    runTest(3, 1, 0.255, 0.015);
}

TEST_CASE("constant_44"){
    runTest(3, 1, 0.5, 0.2);
}

TEST_CASE("constant_45"){
    runTest(3, 2, 0.01, 0.355);
}

TEST_CASE("constant_46"){
    runTest(3, 2, 0.255, 0.015);
}

TEST_CASE("constant_47"){
    runTest(3, 2, 0.5, 0.2);
}

TEST_CASE("constant_48"){
    runTest(2, 1, 0.01, 0.355);
}

TEST_CASE("constant_49"){
    runTest(2, 1, 0.255, 0.015);
}

TEST_CASE("constant_50"){
    runTest(2, 1, 0.5, 0.2);
}

TEST_CASE("constant_51"){
    runTest(2, 2, 0.01, 0.355);
}

TEST_CASE("constant_52"){
    runTest(2, 2, 0.255, 0.015);
}

TEST_CASE("constant_53"){
    runTest(2, 2, 0.5, 0.2);
}

TEST_CASE("constant_54"){
    runTest(3, 1, 0.01, 0.355);
}

TEST_CASE("constant_55"){
    runTest(3, 1, 0.255, 0.015);
}

TEST_CASE("constant_56"){
    runTest(3, 1, 0.5, 0.2);
}

TEST_CASE("constant_57"){
    runTest(3, 2, 0.01, 0.355);
}

TEST_CASE("constant_58"){
    runTest(3, 2, 0.255, 0.015);
}

TEST_CASE("constant_59"){
    runTest(3, 2, 0.5, 0.2);
}

TEST_CASE("constant_60"){
    runTest(2, 1, 0.01, 0.355);
}

TEST_CASE("constant_61"){
    runTest(2, 1, 0.255, 0.015);
}

TEST_CASE("constant_62"){
    runTest(2, 1, 0.5, 0.2);
}

TEST_CASE("constant_63"){
    runTest(2, 2, 0.01, 0.355);
}

TEST_CASE("constant_64"){
    runTest(2, 2, 0.255, 0.015);
}

TEST_CASE("constant_65"){
    runTest(2, 2, 0.5, 0.2);
}

TEST_CASE("constant_66"){
    runTest(3, 1, 0.01, 0.355);
}

TEST_CASE("constant_67"){
    runTest(3, 1, 0.255, 0.015);
}

TEST_CASE("constant_68"){
    runTest(3, 1, 0.5, 0.2);
}

TEST_CASE("constant_69"){
    runTest(3, 2, 0.01, 0.355);
}

TEST_CASE("constant_70"){
    runTest(3, 2, 0.255, 0.015);
}

TEST_CASE("constant_71"){
    runTest(3, 2, 0.5, 0.2);
}




