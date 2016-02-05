/**********************************************************************
 *  Test: 3d inverse fourier transform of constant functions
 *
 *  Analytic source: by hand
 *
 * ********************************************************************/

#include <unittest++/UnitTest++.h>
#include <armadillo>
#include <iostream>

#include "integrator.h"

using namespace std;
using namespace arma;
using namespace edog;

void runTest(int ns, int nt, double dt, double c)
{
    int Ns = pow(2,ns);
    int Nt = pow(2,nt);

    double ds = 0.01;

    Integrator integrator(nt, dt, ns, ds);

    vec k = integrator.spatialFreqVec();
    vec w = integrator.temporalFreqVec();

    cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
    cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
    cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


    //Spatiotemporal signal
    for(int l = 0; l < Nt; l++){
        for(int i = 0; i < Ns; i++){
            for(int j = 0; j < Ns; j++){
                g(i,j,l) = c;
            }
        }
    }

    //fourier signal
    for(int l = 0; l < Nt; l++){
        for(int i = 0; i < Ns; i++){
            for(int j = 0; j < Ns; j++){
                f(i,j,l) = Functions::delta(k[i], 0)
                        * Functions::delta(k[j], 0)
                        * Functions::delta(w[l], 0);
            }
        }
    }

    f *= c*8*PI*PI*PI;
    f /= integrator.spatialFreqResolution()
            * integrator.spatialFreqResolution()
            * integrator.temporalFreqResolution();

    // Backward
    G = integrator.backwardFFT(f);


    // Test
    for(int l = 0; l < Nt; l++){
        for(int i = 0; i < Ns; i++){
            for(int j = 0; j < Ns; j++){
                CHECK_CLOSE(real(g(i,j,l)),
                            real(G(i,j,l)), 1e-9);
            }
        }
    }

}

SUITE(INTEGRATOR){

    TEST(constant_0){
         runTest(2, 1, 0.01, -60.5);
    }

    TEST(constant_1){
         runTest(2, 1, 0.255, -60.5);
    }

    TEST(constant_2){
         runTest(2, 1, 0.5, -60.5);
    }

    TEST(constant_3){
         runTest(2, 2, 0.01, -60.5);
    }

    TEST(constant_4){
         runTest(2, 2, 0.255, -60.5);
    }

    TEST(constant_5){
         runTest(2, 2, 0.5, -60.5);
    }

    TEST(constant_6){
         runTest(3, 1, 0.01, -60.5);
    }

    TEST(constant_7){
         runTest(3, 1, 0.255, -60.5);
    }

    TEST(constant_8){
         runTest(3, 1, 0.5, -60.5);
    }

    TEST(constant_9){
         runTest(3, 2, 0.01, -60.5);
    }

    TEST(constant_10){
         runTest(3, 2, 0.255, -60.5);
    }

    TEST(constant_11){
         runTest(3, 2, 0.5, -60.5);
    }

    TEST(constant_12){
         runTest(2, 1, 0.01, -0.002);
    }

    TEST(constant_13){
         runTest(2, 1, 0.255, -0.002);
    }

    TEST(constant_14){
         runTest(2, 1, 0.5, -0.002);
    }

    TEST(constant_15){
         runTest(2, 2, 0.01, -0.002);
    }

    TEST(constant_16){
         runTest(2, 2, 0.255, -0.002);
    }

    TEST(constant_17){
         runTest(2, 2, 0.5, -0.002);
    }

    TEST(constant_18){
         runTest(3, 1, 0.01, -0.002);
    }

    TEST(constant_19){
         runTest(3, 1, 0.255, -0.002);
    }

    TEST(constant_20){
         runTest(3, 1, 0.5, -0.002);
    }

    TEST(constant_21){
         runTest(3, 2, 0.01, -0.002);
    }

    TEST(constant_22){
         runTest(3, 2, 0.255, -0.002);
    }

    TEST(constant_23){
         runTest(3, 2, 0.5, -0.002);
    }

    TEST(constant_24){
         runTest(2, 1, 0.01, 0.0);
    }

    TEST(constant_25){
         runTest(2, 1, 0.255, 0.0);
    }

    TEST(constant_26){
         runTest(2, 1, 0.5, 0.0);
    }

    TEST(constant_27){
         runTest(2, 2, 0.01, 0.0);
    }

    TEST(constant_28){
         runTest(2, 2, 0.255, 0.0);
    }

    TEST(constant_29){
         runTest(2, 2, 0.5, 0.0);
    }

    TEST(constant_30){
         runTest(3, 1, 0.01, 0.0);
    }

    TEST(constant_31){
         runTest(3, 1, 0.255, 0.0);
    }

    TEST(constant_32){
         runTest(3, 1, 0.5, 0.0);
    }

    TEST(constant_33){
         runTest(3, 2, 0.01, 0.0);
    }

    TEST(constant_34){
         runTest(3, 2, 0.255, 0.0);
    }

    TEST(constant_35){
         runTest(3, 2, 0.5, 0.0);
    }

    TEST(constant_36){
         runTest(2, 1, 0.01, 4.6);
    }

    TEST(constant_37){
         runTest(2, 1, 0.255, 4.6);
    }

    TEST(constant_38){
         runTest(2, 1, 0.5, 4.6);
    }

    TEST(constant_39){
         runTest(2, 2, 0.01, 4.6);
    }

    TEST(constant_40){
         runTest(2, 2, 0.255, 4.6);
    }

    TEST(constant_41){
         runTest(2, 2, 0.5, 4.6);
    }

    TEST(constant_42){
         runTest(3, 1, 0.01, 4.6);
    }

    TEST(constant_43){
         runTest(3, 1, 0.255, 4.6);
    }

    TEST(constant_44){
         runTest(3, 1, 0.5, 4.6);
    }

    TEST(constant_45){
         runTest(3, 2, 0.01, 4.6);
    }

    TEST(constant_46){
         runTest(3, 2, 0.255, 4.6);
    }

    TEST(constant_47){
         runTest(3, 2, 0.5, 4.6);
    }

    TEST(constant_48){
         runTest(2, 1, 0.01, 32.1);
    }

    TEST(constant_49){
         runTest(2, 1, 0.255, 32.1);
    }

    TEST(constant_50){
         runTest(2, 1, 0.5, 32.1);
    }

    TEST(constant_51){
         runTest(2, 2, 0.01, 32.1);
    }

    TEST(constant_52){
         runTest(2, 2, 0.255, 32.1);
    }

    TEST(constant_53){
         runTest(2, 2, 0.5, 32.1);
    }

    TEST(constant_54){
         runTest(3, 1, 0.01, 32.1);
    }

    TEST(constant_55){
         runTest(3, 1, 0.255, 32.1);
    }

    TEST(constant_56){
         runTest(3, 1, 0.5, 32.1);
    }

    TEST(constant_57){
         runTest(3, 2, 0.01, 32.1);
    }

    TEST(constant_58){
         runTest(3, 2, 0.255, 32.1);
    }

    TEST(constant_59){
         runTest(3, 2, 0.5, 32.1);
    }

    TEST(constant_60){
         runTest(2, 1, 0.01, 1000.5);
    }

    TEST(constant_61){
         runTest(2, 1, 0.255, 1000.5);
    }

    TEST(constant_62){
         runTest(2, 1, 0.5, 1000.5);
    }

    TEST(constant_63){
         runTest(2, 2, 0.01, 1000.5);
    }

    TEST(constant_64){
         runTest(2, 2, 0.255, 1000.5);
    }

    TEST(constant_65){
         runTest(2, 2, 0.5, 1000.5);
    }

    TEST(constant_66){
         runTest(3, 1, 0.01, 1000.5);
    }

    TEST(constant_67){
         runTest(3, 1, 0.255, 1000.5);
    }

    TEST(constant_68){
         runTest(3, 1, 0.5, 1000.5);
    }

    TEST(constant_69){
         runTest(3, 2, 0.01, 1000.5);
    }

    TEST(constant_70){
         runTest(3, 2, 0.255, 1000.5);
    }

    TEST(constant_71){
         runTest(3, 2, 0.5, 1000.5);
    }



}



