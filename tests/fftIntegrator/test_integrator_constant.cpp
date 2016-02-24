/**********************************************************************
 *  Test: 3D inverse fourier transform of a constant functions
 *        F(r,t) = F(r) * F(t) = C * C
 *
 *  Analytic source: closed-form experssion
 *
 * ********************************************************************/

#include <unittest++/UnitTest++.h>
#include <lgnSimulator.h>

using namespace lgnSimulator;

void runTest(int ns, int nt, double dt, double ds, double C)
{

    Integrator integrator(nt, dt, ns, ds);
    vec r = integrator.spatialVec();
    vec k = integrator.spatialFreqVec();
    vec t = integrator.timeVec();
    vec w = integrator.temporalFreqVec();

    TemporallyConstant Ft(C, integrator.temporalFreqResolution());
    SpatiallyConstant  Fs(C, integrator.spatialFreqResolution());

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
    cx_cube F_fft = integrator.backwardFFT(G);
    cx_cube diff = F - F_fft;

    cube diff_real = abs(real(diff));
    cube diff_imag = abs(imag(diff));

    //Test
    CHECK_CLOSE(diff_real.max(), 0.0, 1e-9);
    CHECK_CLOSE(diff_imag.max(), 0.0, 1e-9);

}

SUITE(integrator){

    TEST(constant_0){
         runTest(2, 1, 0.01, 0.355, -60.5);
    }

    TEST(constant_1){
         runTest(2, 1, 0.255, 0.015, -60.5);
    }

    TEST(constant_2){
         runTest(2, 1, 0.5, 0.2, -60.5);
    }

    TEST(constant_3){
         runTest(2, 2, 0.01, 0.355, -60.5);
    }

    TEST(constant_4){
         runTest(2, 2, 0.255, 0.015, -60.5);
    }

    TEST(constant_5){
         runTest(2, 2, 0.5, 0.2, -60.5);
    }

    TEST(constant_6){
         runTest(3, 1, 0.01, 0.355, -60.5);
    }

    TEST(constant_7){
         runTest(3, 1, 0.255, 0.015, -60.5);
    }

    TEST(constant_8){
         runTest(3, 1, 0.5, 0.2, -60.5);
    }

    TEST(constant_9){
         runTest(3, 2, 0.01, 0.355, -60.5);
    }

    TEST(constant_10){
         runTest(3, 2, 0.255, 0.015, -60.5);
    }

    TEST(constant_11){
         runTest(3, 2, 0.5, 0.2, -60.5);
    }

    TEST(constant_12){
         runTest(2, 1, 0.01, 0.355, -0.002);
    }

    TEST(constant_13){
         runTest(2, 1, 0.255, 0.015, -0.002);
    }

    TEST(constant_14){
         runTest(2, 1, 0.5, 0.2, -0.002);
    }

    TEST(constant_15){
         runTest(2, 2, 0.01, 0.355, -0.002);
    }

    TEST(constant_16){
         runTest(2, 2, 0.255, 0.015, -0.002);
    }

    TEST(constant_17){
         runTest(2, 2, 0.5, 0.2, -0.002);
    }

    TEST(constant_18){
         runTest(3, 1, 0.01, 0.355, -0.002);
    }

    TEST(constant_19){
         runTest(3, 1, 0.255, 0.015, -0.002);
    }

    TEST(constant_20){
         runTest(3, 1, 0.5, 0.2, -0.002);
    }

    TEST(constant_21){
         runTest(3, 2, 0.01, 0.355, -0.002);
    }

    TEST(constant_22){
         runTest(3, 2, 0.255, 0.015, -0.002);
    }

    TEST(constant_23){
         runTest(3, 2, 0.5, 0.2, -0.002);
    }

    TEST(constant_24){
         runTest(2, 1, 0.01, 0.355, 0.0);
    }

    TEST(constant_25){
         runTest(2, 1, 0.255, 0.015, 0.0);
    }

    TEST(constant_26){
         runTest(2, 1, 0.5, 0.2, 0.0);
    }

    TEST(constant_27){
         runTest(2, 2, 0.01, 0.355, 0.0);
    }

    TEST(constant_28){
         runTest(2, 2, 0.255, 0.015, 0.0);
    }

    TEST(constant_29){
         runTest(2, 2, 0.5, 0.2, 0.0);
    }

    TEST(constant_30){
         runTest(3, 1, 0.01, 0.355, 0.0);
    }

    TEST(constant_31){
         runTest(3, 1, 0.255, 0.015, 0.0);
    }

    TEST(constant_32){
         runTest(3, 1, 0.5, 0.2, 0.0);
    }

    TEST(constant_33){
         runTest(3, 2, 0.01, 0.355, 0.0);
    }

    TEST(constant_34){
         runTest(3, 2, 0.255, 0.015, 0.0);
    }

    TEST(constant_35){
         runTest(3, 2, 0.5, 0.2, 0.0);
    }

    TEST(constant_36){
         runTest(2, 1, 0.01, 0.355, 4.6);
    }

    TEST(constant_37){
         runTest(2, 1, 0.255, 0.015, 4.6);
    }

    TEST(constant_38){
         runTest(2, 1, 0.5, 0.2, 4.6);
    }

    TEST(constant_39){
         runTest(2, 2, 0.01, 0.355, 4.6);
    }

    TEST(constant_40){
         runTest(2, 2, 0.255, 0.015, 4.6);
    }

    TEST(constant_41){
         runTest(2, 2, 0.5, 0.2, 4.6);
    }

    TEST(constant_42){
         runTest(3, 1, 0.01, 0.355, 4.6);
    }

    TEST(constant_43){
         runTest(3, 1, 0.255, 0.015, 4.6);
    }

    TEST(constant_44){
         runTest(3, 1, 0.5, 0.2, 4.6);
    }

    TEST(constant_45){
         runTest(3, 2, 0.01, 0.355, 4.6);
    }

    TEST(constant_46){
         runTest(3, 2, 0.255, 0.015, 4.6);
    }

    TEST(constant_47){
         runTest(3, 2, 0.5, 0.2, 4.6);
    }

    TEST(constant_48){
         runTest(2, 1, 0.01, 0.355, 32.1);
    }

    TEST(constant_49){
         runTest(2, 1, 0.255, 0.015, 32.1);
    }

    TEST(constant_50){
         runTest(2, 1, 0.5, 0.2, 32.1);
    }

    TEST(constant_51){
         runTest(2, 2, 0.01, 0.355, 32.1);
    }

    TEST(constant_52){
         runTest(2, 2, 0.255, 0.015, 32.1);
    }

    TEST(constant_53){
         runTest(2, 2, 0.5, 0.2, 32.1);
    }

    TEST(constant_54){
         runTest(3, 1, 0.01, 0.355, 32.1);
    }

    TEST(constant_55){
         runTest(3, 1, 0.255, 0.015, 32.1);
    }

    TEST(constant_56){
         runTest(3, 1, 0.5, 0.2, 32.1);
    }

    TEST(constant_57){
         runTest(3, 2, 0.01, 0.355, 32.1);
    }

    TEST(constant_58){
         runTest(3, 2, 0.255, 0.015, 32.1);
    }

    TEST(constant_59){
         runTest(3, 2, 0.5, 0.2, 32.1);
    }

    TEST(constant_60){
         runTest(2, 1, 0.01, 0.355, 1000.5);
    }

    TEST(constant_61){
         runTest(2, 1, 0.255, 0.015, 1000.5);
    }

    TEST(constant_62){
         runTest(2, 1, 0.5, 0.2, 1000.5);
    }

    TEST(constant_63){
         runTest(2, 2, 0.01, 0.355, 1000.5);
    }

    TEST(constant_64){
         runTest(2, 2, 0.255, 0.015, 1000.5);
    }

    TEST(constant_65){
         runTest(2, 2, 0.5, 0.2, 1000.5);
    }

    TEST(constant_66){
         runTest(3, 1, 0.01, 0.355, 1000.5);
    }

    TEST(constant_67){
         runTest(3, 1, 0.255, 0.015, 1000.5);
    }

    TEST(constant_68){
         runTest(3, 1, 0.5, 0.2, 1000.5);
    }

    TEST(constant_69){
         runTest(3, 2, 0.01, 0.355, 1000.5);
    }

    TEST(constant_70){
         runTest(3, 2, 0.255, 0.015, 1000.5);
    }

    TEST(constant_71){
         runTest(3, 2, 0.5, 0.2, 1000.5);
    }



}



