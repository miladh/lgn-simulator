/**********************************************************************
 *  Test: 3d inverse fourier transform of cosine waves
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
using namespace lgnSimulator;

void runTest(int ns, int nt, double dt, int wdId, int kxId, int kyId)
{
    int Ns = pow(2,ns);
    int Nt = pow(2,nt);

    double ds = 0.01;

    Integrator integrator(nt, dt, ns, ds);

    vec s = integrator.spatialVec();
    vec k = integrator.spatialFreqVec();
    vec t = integrator.timeVec();
    vec w = integrator.temporalFreqVec();

    cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
    cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
    cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);

    double wd = w(wdId);
    double kx = k(kxId);
    double ky = k(kyId);

    //Spatiotemporal signal
    for(int l = 0; l < Nt; l++){
      for(int i = 0; i < Ns; i++){
        for(int j = 0; j < Ns; j++){
          g(i,j,l) = cos(kx*s[i] + ky*s[j] - wd * t[l]);
        }
      }
    }

    //fourier signal
    for(int l = 0; l < Nt; l++){
      for(int i = 0; i < Ns; i++){
        for(int j = 0; j < Ns; j++){
          f(i,j,l) = Special::delta(k[i], kx)
                   * Special::delta(k[j], ky)
                   * Special::delta(w[l], -wd);
        }
      }
    }

    f *= 8*core::pi*core::pi*core::pi;
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

SUITE(integrator){

    TEST(cosine_0){
         runTest(2, 2, 0.01, 0, 1, 1);
    }

    TEST(cosine_1){
         runTest(2, 2, 0.01, 0, 1, 3);
    }

    TEST(cosine_2){
         runTest(2, 2, 0.01, 0, 3, 1);
    }

    TEST(cosine_3){
         runTest(2, 2, 0.01, 0, 3, 3);
    }

    TEST(cosine_4){
         runTest(2, 2, 0.01, 1, 1, 1);
    }

    TEST(cosine_5){
         runTest(2, 2, 0.01, 1, 1, 3);
    }

    TEST(cosine_6){
         runTest(2, 2, 0.01, 1, 3, 1);
    }

    TEST(cosine_7){
         runTest(2, 2, 0.01, 1, 3, 3);
    }

    TEST(cosine_8){
         runTest(2, 2, 0.01, 3, 1, 1);
    }

    TEST(cosine_9){
         runTest(2, 2, 0.01, 3, 1, 3);
    }

    TEST(cosine_10){
         runTest(2, 2, 0.01, 3, 3, 1);
    }

    TEST(cosine_11){
         runTest(2, 2, 0.01, 3, 3, 3);
    }

    TEST(cosine_12){
         runTest(2, 2, 0.5, 0, 1, 1);
    }

    TEST(cosine_13){
         runTest(2, 2, 0.5, 0, 1, 3);
    }

    TEST(cosine_14){
         runTest(2, 2, 0.5, 0, 3, 1);
    }

    TEST(cosine_15){
         runTest(2, 2, 0.5, 0, 3, 3);
    }

    TEST(cosine_16){
         runTest(2, 2, 0.5, 1, 1, 1);
    }

    TEST(cosine_17){
         runTest(2, 2, 0.5, 1, 1, 3);
    }

    TEST(cosine_18){
         runTest(2, 2, 0.5, 1, 3, 1);
    }

    TEST(cosine_19){
         runTest(2, 2, 0.5, 1, 3, 3);
    }

    TEST(cosine_20){
         runTest(2, 2, 0.5, 3, 1, 1);
    }

    TEST(cosine_21){
         runTest(2, 2, 0.5, 3, 1, 3);
    }

    TEST(cosine_22){
         runTest(2, 2, 0.5, 3, 3, 1);
    }

    TEST(cosine_23){
         runTest(2, 2, 0.5, 3, 3, 3);
    }

    TEST(cosine_24){
         runTest(2, 3, 0.01, 0, 1, 1);
    }

    TEST(cosine_25){
         runTest(2, 3, 0.01, 0, 1, 3);
    }

    TEST(cosine_26){
         runTest(2, 3, 0.01, 0, 3, 1);
    }

    TEST(cosine_27){
         runTest(2, 3, 0.01, 0, 3, 3);
    }

    TEST(cosine_28){
         runTest(2, 3, 0.01, 1, 1, 1);
    }

    TEST(cosine_29){
         runTest(2, 3, 0.01, 1, 1, 3);
    }

    TEST(cosine_30){
         runTest(2, 3, 0.01, 1, 3, 1);
    }

    TEST(cosine_31){
         runTest(2, 3, 0.01, 1, 3, 3);
    }

    TEST(cosine_32){
         runTest(2, 3, 0.01, 2, 1, 1);
    }

    TEST(cosine_33){
         runTest(2, 3, 0.01, 2, 1, 3);
    }

    TEST(cosine_34){
         runTest(2, 3, 0.01, 2, 3, 1);
    }

    TEST(cosine_35){
         runTest(2, 3, 0.01, 2, 3, 3);
    }

    TEST(cosine_36){
         runTest(2, 3, 0.01, 3, 1, 1);
    }

    TEST(cosine_37){
         runTest(2, 3, 0.01, 3, 1, 3);
    }

    TEST(cosine_38){
         runTest(2, 3, 0.01, 3, 3, 1);
    }

    TEST(cosine_39){
         runTest(2, 3, 0.01, 3, 3, 3);
    }

    TEST(cosine_40){
         runTest(2, 3, 0.01, 5, 1, 1);
    }

    TEST(cosine_41){
         runTest(2, 3, 0.01, 5, 1, 3);
    }

    TEST(cosine_42){
         runTest(2, 3, 0.01, 5, 3, 1);
    }

    TEST(cosine_43){
         runTest(2, 3, 0.01, 5, 3, 3);
    }

    TEST(cosine_44){
         runTest(2, 3, 0.01, 6, 1, 1);
    }

    TEST(cosine_45){
         runTest(2, 3, 0.01, 6, 1, 3);
    }

    TEST(cosine_46){
         runTest(2, 3, 0.01, 6, 3, 1);
    }

    TEST(cosine_47){
         runTest(2, 3, 0.01, 6, 3, 3);
    }

    TEST(cosine_48){
         runTest(2, 3, 0.01, 7, 1, 1);
    }

    TEST(cosine_49){
         runTest(2, 3, 0.01, 7, 1, 3);
    }

    TEST(cosine_50){
         runTest(2, 3, 0.01, 7, 3, 1);
    }

    TEST(cosine_51){
         runTest(2, 3, 0.01, 7, 3, 3);
    }

    TEST(cosine_52){
         runTest(2, 3, 0.5, 0, 1, 1);
    }

    TEST(cosine_53){
         runTest(2, 3, 0.5, 0, 1, 3);
    }

    TEST(cosine_54){
         runTest(2, 3, 0.5, 0, 3, 1);
    }

    TEST(cosine_55){
         runTest(2, 3, 0.5, 0, 3, 3);
    }

    TEST(cosine_56){
         runTest(2, 3, 0.5, 1, 1, 1);
    }

    TEST(cosine_57){
         runTest(2, 3, 0.5, 1, 1, 3);
    }

    TEST(cosine_58){
         runTest(2, 3, 0.5, 1, 3, 1);
    }

    TEST(cosine_59){
         runTest(2, 3, 0.5, 1, 3, 3);
    }

    TEST(cosine_60){
         runTest(2, 3, 0.5, 2, 1, 1);
    }

    TEST(cosine_61){
         runTest(2, 3, 0.5, 2, 1, 3);
    }

    TEST(cosine_62){
         runTest(2, 3, 0.5, 2, 3, 1);
    }

    TEST(cosine_63){
         runTest(2, 3, 0.5, 2, 3, 3);
    }

    TEST(cosine_64){
         runTest(2, 3, 0.5, 3, 1, 1);
    }

    TEST(cosine_65){
         runTest(2, 3, 0.5, 3, 1, 3);
    }

    TEST(cosine_66){
         runTest(2, 3, 0.5, 3, 3, 1);
    }

    TEST(cosine_67){
         runTest(2, 3, 0.5, 3, 3, 3);
    }

    TEST(cosine_68){
         runTest(2, 3, 0.5, 5, 1, 1);
    }

    TEST(cosine_69){
         runTest(2, 3, 0.5, 5, 1, 3);
    }

    TEST(cosine_70){
         runTest(2, 3, 0.5, 5, 3, 1);
    }

    TEST(cosine_71){
         runTest(2, 3, 0.5, 5, 3, 3);
    }

    TEST(cosine_72){
         runTest(2, 3, 0.5, 6, 1, 1);
    }

    TEST(cosine_73){
         runTest(2, 3, 0.5, 6, 1, 3);
    }

    TEST(cosine_74){
         runTest(2, 3, 0.5, 6, 3, 1);
    }

    TEST(cosine_75){
         runTest(2, 3, 0.5, 6, 3, 3);
    }

    TEST(cosine_76){
         runTest(2, 3, 0.5, 7, 1, 1);
    }

    TEST(cosine_77){
         runTest(2, 3, 0.5, 7, 1, 3);
    }

    TEST(cosine_78){
         runTest(2, 3, 0.5, 7, 3, 1);
    }

    TEST(cosine_79){
         runTest(2, 3, 0.5, 7, 3, 3);
    }

    TEST(cosine_80){
         runTest(3, 2, 0.01, 0, 1, 1);
    }

    TEST(cosine_81){
         runTest(3, 2, 0.01, 0, 1, 2);
    }

    TEST(cosine_82){
         runTest(3, 2, 0.01, 0, 1, 3);
    }

    TEST(cosine_83){
         runTest(3, 2, 0.01, 0, 1, 5);
    }

    TEST(cosine_84){
         runTest(3, 2, 0.01, 0, 1, 6);
    }

    TEST(cosine_85){
         runTest(3, 2, 0.01, 0, 1, 7);
    }

    TEST(cosine_86){
         runTest(3, 2, 0.01, 0, 2, 1);
    }

    TEST(cosine_87){
         runTest(3, 2, 0.01, 0, 2, 2);
    }

    TEST(cosine_88){
         runTest(3, 2, 0.01, 0, 2, 3);
    }

    TEST(cosine_89){
         runTest(3, 2, 0.01, 0, 2, 5);
    }

    TEST(cosine_90){
         runTest(3, 2, 0.01, 0, 2, 6);
    }

    TEST(cosine_91){
         runTest(3, 2, 0.01, 0, 2, 7);
    }

    TEST(cosine_92){
         runTest(3, 2, 0.01, 0, 3, 1);
    }

    TEST(cosine_93){
         runTest(3, 2, 0.01, 0, 3, 2);
    }

    TEST(cosine_94){
         runTest(3, 2, 0.01, 0, 3, 3);
    }

    TEST(cosine_95){
         runTest(3, 2, 0.01, 0, 3, 5);
    }

    TEST(cosine_96){
         runTest(3, 2, 0.01, 0, 3, 6);
    }

    TEST(cosine_97){
         runTest(3, 2, 0.01, 0, 3, 7);
    }

    TEST(cosine_98){
         runTest(3, 2, 0.01, 0, 5, 1);
    }

    TEST(cosine_99){
         runTest(3, 2, 0.01, 0, 5, 2);
    }

    TEST(cosine_100){
         runTest(3, 2, 0.01, 0, 5, 3);
    }

    TEST(cosine_101){
         runTest(3, 2, 0.01, 0, 5, 5);
    }

    TEST(cosine_102){
         runTest(3, 2, 0.01, 0, 5, 6);
    }

    TEST(cosine_103){
         runTest(3, 2, 0.01, 0, 5, 7);
    }

    TEST(cosine_104){
         runTest(3, 2, 0.01, 0, 6, 1);
    }

    TEST(cosine_105){
         runTest(3, 2, 0.01, 0, 6, 2);
    }

    TEST(cosine_106){
         runTest(3, 2, 0.01, 0, 6, 3);
    }

    TEST(cosine_107){
         runTest(3, 2, 0.01, 0, 6, 5);
    }

    TEST(cosine_108){
         runTest(3, 2, 0.01, 0, 6, 6);
    }

    TEST(cosine_109){
         runTest(3, 2, 0.01, 0, 6, 7);
    }

    TEST(cosine_110){
         runTest(3, 2, 0.01, 0, 7, 1);
    }

    TEST(cosine_111){
         runTest(3, 2, 0.01, 0, 7, 2);
    }

    TEST(cosine_112){
         runTest(3, 2, 0.01, 0, 7, 3);
    }

    TEST(cosine_113){
         runTest(3, 2, 0.01, 0, 7, 5);
    }

    TEST(cosine_114){
         runTest(3, 2, 0.01, 0, 7, 6);
    }

    TEST(cosine_115){
         runTest(3, 2, 0.01, 0, 7, 7);
    }

    TEST(cosine_116){
         runTest(3, 2, 0.01, 1, 1, 1);
    }

    TEST(cosine_117){
         runTest(3, 2, 0.01, 1, 1, 2);
    }

    TEST(cosine_118){
         runTest(3, 2, 0.01, 1, 1, 3);
    }

    TEST(cosine_119){
         runTest(3, 2, 0.01, 1, 1, 5);
    }

    TEST(cosine_120){
         runTest(3, 2, 0.01, 1, 1, 6);
    }

    TEST(cosine_121){
         runTest(3, 2, 0.01, 1, 1, 7);
    }

    TEST(cosine_122){
         runTest(3, 2, 0.01, 1, 2, 1);
    }

    TEST(cosine_123){
         runTest(3, 2, 0.01, 1, 2, 2);
    }

    TEST(cosine_124){
         runTest(3, 2, 0.01, 1, 2, 3);
    }

    TEST(cosine_125){
         runTest(3, 2, 0.01, 1, 2, 5);
    }

    TEST(cosine_126){
         runTest(3, 2, 0.01, 1, 2, 6);
    }

    TEST(cosine_127){
         runTest(3, 2, 0.01, 1, 2, 7);
    }

    TEST(cosine_128){
         runTest(3, 2, 0.01, 1, 3, 1);
    }

    TEST(cosine_129){
         runTest(3, 2, 0.01, 1, 3, 2);
    }

    TEST(cosine_130){
         runTest(3, 2, 0.01, 1, 3, 3);
    }

    TEST(cosine_131){
         runTest(3, 2, 0.01, 1, 3, 5);
    }

    TEST(cosine_132){
         runTest(3, 2, 0.01, 1, 3, 6);
    }

    TEST(cosine_133){
         runTest(3, 2, 0.01, 1, 3, 7);
    }

    TEST(cosine_134){
         runTest(3, 2, 0.01, 1, 5, 1);
    }

    TEST(cosine_135){
         runTest(3, 2, 0.01, 1, 5, 2);
    }

    TEST(cosine_136){
         runTest(3, 2, 0.01, 1, 5, 3);
    }

    TEST(cosine_137){
         runTest(3, 2, 0.01, 1, 5, 5);
    }

    TEST(cosine_138){
         runTest(3, 2, 0.01, 1, 5, 6);
    }

    TEST(cosine_139){
         runTest(3, 2, 0.01, 1, 5, 7);
    }

    TEST(cosine_140){
         runTest(3, 2, 0.01, 1, 6, 1);
    }

    TEST(cosine_141){
         runTest(3, 2, 0.01, 1, 6, 2);
    }

    TEST(cosine_142){
         runTest(3, 2, 0.01, 1, 6, 3);
    }

    TEST(cosine_143){
         runTest(3, 2, 0.01, 1, 6, 5);
    }

    TEST(cosine_144){
         runTest(3, 2, 0.01, 1, 6, 6);
    }

    TEST(cosine_145){
         runTest(3, 2, 0.01, 1, 6, 7);
    }

    TEST(cosine_146){
         runTest(3, 2, 0.01, 1, 7, 1);
    }

    TEST(cosine_147){
         runTest(3, 2, 0.01, 1, 7, 2);
    }

    TEST(cosine_148){
         runTest(3, 2, 0.01, 1, 7, 3);
    }

    TEST(cosine_149){
         runTest(3, 2, 0.01, 1, 7, 5);
    }

    TEST(cosine_150){
         runTest(3, 2, 0.01, 1, 7, 6);
    }

    TEST(cosine_151){
         runTest(3, 2, 0.01, 1, 7, 7);
    }

    TEST(cosine_152){
         runTest(3, 2, 0.01, 3, 1, 1);
    }

    TEST(cosine_153){
         runTest(3, 2, 0.01, 3, 1, 2);
    }

    TEST(cosine_154){
         runTest(3, 2, 0.01, 3, 1, 3);
    }

    TEST(cosine_155){
         runTest(3, 2, 0.01, 3, 1, 5);
    }

    TEST(cosine_156){
         runTest(3, 2, 0.01, 3, 1, 6);
    }

    TEST(cosine_157){
         runTest(3, 2, 0.01, 3, 1, 7);
    }

    TEST(cosine_158){
         runTest(3, 2, 0.01, 3, 2, 1);
    }

    TEST(cosine_159){
         runTest(3, 2, 0.01, 3, 2, 2);
    }

    TEST(cosine_160){
         runTest(3, 2, 0.01, 3, 2, 3);
    }

    TEST(cosine_161){
         runTest(3, 2, 0.01, 3, 2, 5);
    }

    TEST(cosine_162){
         runTest(3, 2, 0.01, 3, 2, 6);
    }

    TEST(cosine_163){
         runTest(3, 2, 0.01, 3, 2, 7);
    }

    TEST(cosine_164){
         runTest(3, 2, 0.01, 3, 3, 1);
    }

    TEST(cosine_165){
         runTest(3, 2, 0.01, 3, 3, 2);
    }

    TEST(cosine_166){
         runTest(3, 2, 0.01, 3, 3, 3);
    }

    TEST(cosine_167){
         runTest(3, 2, 0.01, 3, 3, 5);
    }

    TEST(cosine_168){
         runTest(3, 2, 0.01, 3, 3, 6);
    }

    TEST(cosine_169){
         runTest(3, 2, 0.01, 3, 3, 7);
    }

    TEST(cosine_170){
         runTest(3, 2, 0.01, 3, 5, 1);
    }

    TEST(cosine_171){
         runTest(3, 2, 0.01, 3, 5, 2);
    }

    TEST(cosine_172){
         runTest(3, 2, 0.01, 3, 5, 3);
    }

    TEST(cosine_173){
         runTest(3, 2, 0.01, 3, 5, 5);
    }

    TEST(cosine_174){
         runTest(3, 2, 0.01, 3, 5, 6);
    }

    TEST(cosine_175){
         runTest(3, 2, 0.01, 3, 5, 7);
    }

    TEST(cosine_176){
         runTest(3, 2, 0.01, 3, 6, 1);
    }

    TEST(cosine_177){
         runTest(3, 2, 0.01, 3, 6, 2);
    }

    TEST(cosine_178){
         runTest(3, 2, 0.01, 3, 6, 3);
    }

    TEST(cosine_179){
         runTest(3, 2, 0.01, 3, 6, 5);
    }

    TEST(cosine_180){
         runTest(3, 2, 0.01, 3, 6, 6);
    }

    TEST(cosine_181){
         runTest(3, 2, 0.01, 3, 6, 7);
    }

    TEST(cosine_182){
         runTest(3, 2, 0.01, 3, 7, 1);
    }

    TEST(cosine_183){
         runTest(3, 2, 0.01, 3, 7, 2);
    }

    TEST(cosine_184){
         runTest(3, 2, 0.01, 3, 7, 3);
    }

    TEST(cosine_185){
         runTest(3, 2, 0.01, 3, 7, 5);
    }

    TEST(cosine_186){
         runTest(3, 2, 0.01, 3, 7, 6);
    }

    TEST(cosine_187){
         runTest(3, 2, 0.01, 3, 7, 7);
    }

    TEST(cosine_188){
         runTest(3, 2, 0.5, 0, 1, 1);
    }

    TEST(cosine_189){
         runTest(3, 2, 0.5, 0, 1, 2);
    }

    TEST(cosine_190){
         runTest(3, 2, 0.5, 0, 1, 3);
    }

    TEST(cosine_191){
         runTest(3, 2, 0.5, 0, 1, 5);
    }

    TEST(cosine_192){
         runTest(3, 2, 0.5, 0, 1, 6);
    }

    TEST(cosine_193){
         runTest(3, 2, 0.5, 0, 1, 7);
    }

    TEST(cosine_194){
         runTest(3, 2, 0.5, 0, 2, 1);
    }

    TEST(cosine_195){
         runTest(3, 2, 0.5, 0, 2, 2);
    }

    TEST(cosine_196){
         runTest(3, 2, 0.5, 0, 2, 3);
    }

    TEST(cosine_197){
         runTest(3, 2, 0.5, 0, 2, 5);
    }

    TEST(cosine_198){
         runTest(3, 2, 0.5, 0, 2, 6);
    }

    TEST(cosine_199){
         runTest(3, 2, 0.5, 0, 2, 7);
    }

    TEST(cosine_200){
         runTest(3, 2, 0.5, 0, 3, 1);
    }

    TEST(cosine_201){
         runTest(3, 2, 0.5, 0, 3, 2);
    }

    TEST(cosine_202){
         runTest(3, 2, 0.5, 0, 3, 3);
    }

    TEST(cosine_203){
         runTest(3, 2, 0.5, 0, 3, 5);
    }

    TEST(cosine_204){
         runTest(3, 2, 0.5, 0, 3, 6);
    }

    TEST(cosine_205){
         runTest(3, 2, 0.5, 0, 3, 7);
    }

    TEST(cosine_206){
         runTest(3, 2, 0.5, 0, 5, 1);
    }

    TEST(cosine_207){
         runTest(3, 2, 0.5, 0, 5, 2);
    }

    TEST(cosine_208){
         runTest(3, 2, 0.5, 0, 5, 3);
    }

    TEST(cosine_209){
         runTest(3, 2, 0.5, 0, 5, 5);
    }

    TEST(cosine_210){
         runTest(3, 2, 0.5, 0, 5, 6);
    }

    TEST(cosine_211){
         runTest(3, 2, 0.5, 0, 5, 7);
    }

    TEST(cosine_212){
         runTest(3, 2, 0.5, 0, 6, 1);
    }

    TEST(cosine_213){
         runTest(3, 2, 0.5, 0, 6, 2);
    }

    TEST(cosine_214){
         runTest(3, 2, 0.5, 0, 6, 3);
    }

    TEST(cosine_215){
         runTest(3, 2, 0.5, 0, 6, 5);
    }

    TEST(cosine_216){
         runTest(3, 2, 0.5, 0, 6, 6);
    }

    TEST(cosine_217){
         runTest(3, 2, 0.5, 0, 6, 7);
    }

    TEST(cosine_218){
         runTest(3, 2, 0.5, 0, 7, 1);
    }

    TEST(cosine_219){
         runTest(3, 2, 0.5, 0, 7, 2);
    }

    TEST(cosine_220){
         runTest(3, 2, 0.5, 0, 7, 3);
    }

    TEST(cosine_221){
         runTest(3, 2, 0.5, 0, 7, 5);
    }

    TEST(cosine_222){
         runTest(3, 2, 0.5, 0, 7, 6);
    }

    TEST(cosine_223){
         runTest(3, 2, 0.5, 0, 7, 7);
    }

    TEST(cosine_224){
         runTest(3, 2, 0.5, 1, 1, 1);
    }

    TEST(cosine_225){
         runTest(3, 2, 0.5, 1, 1, 2);
    }

    TEST(cosine_226){
         runTest(3, 2, 0.5, 1, 1, 3);
    }

    TEST(cosine_227){
         runTest(3, 2, 0.5, 1, 1, 5);
    }

    TEST(cosine_228){
         runTest(3, 2, 0.5, 1, 1, 6);
    }

    TEST(cosine_229){
         runTest(3, 2, 0.5, 1, 1, 7);
    }

    TEST(cosine_230){
         runTest(3, 2, 0.5, 1, 2, 1);
    }

    TEST(cosine_231){
         runTest(3, 2, 0.5, 1, 2, 2);
    }

    TEST(cosine_232){
         runTest(3, 2, 0.5, 1, 2, 3);
    }

    TEST(cosine_233){
         runTest(3, 2, 0.5, 1, 2, 5);
    }

    TEST(cosine_234){
         runTest(3, 2, 0.5, 1, 2, 6);
    }

    TEST(cosine_235){
         runTest(3, 2, 0.5, 1, 2, 7);
    }

    TEST(cosine_236){
         runTest(3, 2, 0.5, 1, 3, 1);
    }

    TEST(cosine_237){
         runTest(3, 2, 0.5, 1, 3, 2);
    }

    TEST(cosine_238){
         runTest(3, 2, 0.5, 1, 3, 3);
    }

    TEST(cosine_239){
         runTest(3, 2, 0.5, 1, 3, 5);
    }

    TEST(cosine_240){
         runTest(3, 2, 0.5, 1, 3, 6);
    }

    TEST(cosine_241){
         runTest(3, 2, 0.5, 1, 3, 7);
    }

    TEST(cosine_242){
         runTest(3, 2, 0.5, 1, 5, 1);
    }

    TEST(cosine_243){
         runTest(3, 2, 0.5, 1, 5, 2);
    }

    TEST(cosine_244){
         runTest(3, 2, 0.5, 1, 5, 3);
    }

    TEST(cosine_245){
         runTest(3, 2, 0.5, 1, 5, 5);
    }

    TEST(cosine_246){
         runTest(3, 2, 0.5, 1, 5, 6);
    }

    TEST(cosine_247){
         runTest(3, 2, 0.5, 1, 5, 7);
    }

    TEST(cosine_248){
         runTest(3, 2, 0.5, 1, 6, 1);
    }

    TEST(cosine_249){
         runTest(3, 2, 0.5, 1, 6, 2);
    }

    TEST(cosine_250){
         runTest(3, 2, 0.5, 1, 6, 3);
    }

    TEST(cosine_251){
         runTest(3, 2, 0.5, 1, 6, 5);
    }

    TEST(cosine_252){
         runTest(3, 2, 0.5, 1, 6, 6);
    }

    TEST(cosine_253){
         runTest(3, 2, 0.5, 1, 6, 7);
    }

    TEST(cosine_254){
         runTest(3, 2, 0.5, 1, 7, 1);
    }

    TEST(cosine_255){
         runTest(3, 2, 0.5, 1, 7, 2);
    }

    TEST(cosine_256){
         runTest(3, 2, 0.5, 1, 7, 3);
    }

    TEST(cosine_257){
         runTest(3, 2, 0.5, 1, 7, 5);
    }

    TEST(cosine_258){
         runTest(3, 2, 0.5, 1, 7, 6);
    }

    TEST(cosine_259){
         runTest(3, 2, 0.5, 1, 7, 7);
    }

    TEST(cosine_260){
         runTest(3, 2, 0.5, 3, 1, 1);
    }

    TEST(cosine_261){
         runTest(3, 2, 0.5, 3, 1, 2);
    }

    TEST(cosine_262){
         runTest(3, 2, 0.5, 3, 1, 3);
    }

    TEST(cosine_263){
         runTest(3, 2, 0.5, 3, 1, 5);
    }

    TEST(cosine_264){
         runTest(3, 2, 0.5, 3, 1, 6);
    }

    TEST(cosine_265){
         runTest(3, 2, 0.5, 3, 1, 7);
    }

    TEST(cosine_266){
         runTest(3, 2, 0.5, 3, 2, 1);
    }

    TEST(cosine_267){
         runTest(3, 2, 0.5, 3, 2, 2);
    }

    TEST(cosine_268){
         runTest(3, 2, 0.5, 3, 2, 3);
    }

    TEST(cosine_269){
         runTest(3, 2, 0.5, 3, 2, 5);
    }

    TEST(cosine_270){
         runTest(3, 2, 0.5, 3, 2, 6);
    }

    TEST(cosine_271){
         runTest(3, 2, 0.5, 3, 2, 7);
    }

    TEST(cosine_272){
         runTest(3, 2, 0.5, 3, 3, 1);
    }

    TEST(cosine_273){
         runTest(3, 2, 0.5, 3, 3, 2);
    }

    TEST(cosine_274){
         runTest(3, 2, 0.5, 3, 3, 3);
    }

    TEST(cosine_275){
         runTest(3, 2, 0.5, 3, 3, 5);
    }

    TEST(cosine_276){
         runTest(3, 2, 0.5, 3, 3, 6);
    }

    TEST(cosine_277){
         runTest(3, 2, 0.5, 3, 3, 7);
    }

    TEST(cosine_278){
         runTest(3, 2, 0.5, 3, 5, 1);
    }

    TEST(cosine_279){
         runTest(3, 2, 0.5, 3, 5, 2);
    }

    TEST(cosine_280){
         runTest(3, 2, 0.5, 3, 5, 3);
    }

    TEST(cosine_281){
         runTest(3, 2, 0.5, 3, 5, 5);
    }

    TEST(cosine_282){
         runTest(3, 2, 0.5, 3, 5, 6);
    }

    TEST(cosine_283){
         runTest(3, 2, 0.5, 3, 5, 7);
    }

    TEST(cosine_284){
         runTest(3, 2, 0.5, 3, 6, 1);
    }

    TEST(cosine_285){
         runTest(3, 2, 0.5, 3, 6, 2);
    }

    TEST(cosine_286){
         runTest(3, 2, 0.5, 3, 6, 3);
    }

    TEST(cosine_287){
         runTest(3, 2, 0.5, 3, 6, 5);
    }

    TEST(cosine_288){
         runTest(3, 2, 0.5, 3, 6, 6);
    }

    TEST(cosine_289){
         runTest(3, 2, 0.5, 3, 6, 7);
    }

    TEST(cosine_290){
         runTest(3, 2, 0.5, 3, 7, 1);
    }

    TEST(cosine_291){
         runTest(3, 2, 0.5, 3, 7, 2);
    }

    TEST(cosine_292){
         runTest(3, 2, 0.5, 3, 7, 3);
    }

    TEST(cosine_293){
         runTest(3, 2, 0.5, 3, 7, 5);
    }

    TEST(cosine_294){
         runTest(3, 2, 0.5, 3, 7, 6);
    }

    TEST(cosine_295){
         runTest(3, 2, 0.5, 3, 7, 7);
    }

    TEST(cosine_296){
         runTest(3, 3, 0.01, 0, 1, 1);
    }

    TEST(cosine_297){
         runTest(3, 3, 0.01, 0, 1, 2);
    }

    TEST(cosine_298){
         runTest(3, 3, 0.01, 0, 1, 3);
    }

    TEST(cosine_299){
         runTest(3, 3, 0.01, 0, 1, 5);
    }

    TEST(cosine_300){
         runTest(3, 3, 0.01, 0, 1, 6);
    }

    TEST(cosine_301){
         runTest(3, 3, 0.01, 0, 1, 7);
    }

    TEST(cosine_302){
         runTest(3, 3, 0.01, 0, 2, 1);
    }

    TEST(cosine_303){
         runTest(3, 3, 0.01, 0, 2, 2);
    }

    TEST(cosine_304){
         runTest(3, 3, 0.01, 0, 2, 3);
    }

    TEST(cosine_305){
         runTest(3, 3, 0.01, 0, 2, 5);
    }

    TEST(cosine_306){
         runTest(3, 3, 0.01, 0, 2, 6);
    }

    TEST(cosine_307){
         runTest(3, 3, 0.01, 0, 2, 7);
    }

    TEST(cosine_308){
         runTest(3, 3, 0.01, 0, 3, 1);
    }

    TEST(cosine_309){
         runTest(3, 3, 0.01, 0, 3, 2);
    }

    TEST(cosine_310){
         runTest(3, 3, 0.01, 0, 3, 3);
    }

    TEST(cosine_311){
         runTest(3, 3, 0.01, 0, 3, 5);
    }

    TEST(cosine_312){
         runTest(3, 3, 0.01, 0, 3, 6);
    }

    TEST(cosine_313){
         runTest(3, 3, 0.01, 0, 3, 7);
    }

    TEST(cosine_314){
         runTest(3, 3, 0.01, 0, 5, 1);
    }

    TEST(cosine_315){
         runTest(3, 3, 0.01, 0, 5, 2);
    }

    TEST(cosine_316){
         runTest(3, 3, 0.01, 0, 5, 3);
    }

    TEST(cosine_317){
         runTest(3, 3, 0.01, 0, 5, 5);
    }

    TEST(cosine_318){
         runTest(3, 3, 0.01, 0, 5, 6);
    }

    TEST(cosine_319){
         runTest(3, 3, 0.01, 0, 5, 7);
    }

    TEST(cosine_320){
         runTest(3, 3, 0.01, 0, 6, 1);
    }

    TEST(cosine_321){
         runTest(3, 3, 0.01, 0, 6, 2);
    }

    TEST(cosine_322){
         runTest(3, 3, 0.01, 0, 6, 3);
    }

    TEST(cosine_323){
         runTest(3, 3, 0.01, 0, 6, 5);
    }

    TEST(cosine_324){
         runTest(3, 3, 0.01, 0, 6, 6);
    }

    TEST(cosine_325){
         runTest(3, 3, 0.01, 0, 6, 7);
    }

    TEST(cosine_326){
         runTest(3, 3, 0.01, 0, 7, 1);
    }

    TEST(cosine_327){
         runTest(3, 3, 0.01, 0, 7, 2);
    }

    TEST(cosine_328){
         runTest(3, 3, 0.01, 0, 7, 3);
    }

    TEST(cosine_329){
         runTest(3, 3, 0.01, 0, 7, 5);
    }

    TEST(cosine_330){
         runTest(3, 3, 0.01, 0, 7, 6);
    }

    TEST(cosine_331){
         runTest(3, 3, 0.01, 0, 7, 7);
    }

    TEST(cosine_332){
         runTest(3, 3, 0.01, 1, 1, 1);
    }

    TEST(cosine_333){
         runTest(3, 3, 0.01, 1, 1, 2);
    }

    TEST(cosine_334){
         runTest(3, 3, 0.01, 1, 1, 3);
    }

    TEST(cosine_335){
         runTest(3, 3, 0.01, 1, 1, 5);
    }

    TEST(cosine_336){
         runTest(3, 3, 0.01, 1, 1, 6);
    }

    TEST(cosine_337){
         runTest(3, 3, 0.01, 1, 1, 7);
    }

    TEST(cosine_338){
         runTest(3, 3, 0.01, 1, 2, 1);
    }

    TEST(cosine_339){
         runTest(3, 3, 0.01, 1, 2, 2);
    }

    TEST(cosine_340){
         runTest(3, 3, 0.01, 1, 2, 3);
    }

    TEST(cosine_341){
         runTest(3, 3, 0.01, 1, 2, 5);
    }

    TEST(cosine_342){
         runTest(3, 3, 0.01, 1, 2, 6);
    }

    TEST(cosine_343){
         runTest(3, 3, 0.01, 1, 2, 7);
    }

    TEST(cosine_344){
         runTest(3, 3, 0.01, 1, 3, 1);
    }

    TEST(cosine_345){
         runTest(3, 3, 0.01, 1, 3, 2);
    }

    TEST(cosine_346){
         runTest(3, 3, 0.01, 1, 3, 3);
    }

    TEST(cosine_347){
         runTest(3, 3, 0.01, 1, 3, 5);
    }

    TEST(cosine_348){
         runTest(3, 3, 0.01, 1, 3, 6);
    }

    TEST(cosine_349){
         runTest(3, 3, 0.01, 1, 3, 7);
    }

    TEST(cosine_350){
         runTest(3, 3, 0.01, 1, 5, 1);
    }

    TEST(cosine_351){
         runTest(3, 3, 0.01, 1, 5, 2);
    }

    TEST(cosine_352){
         runTest(3, 3, 0.01, 1, 5, 3);
    }

    TEST(cosine_353){
         runTest(3, 3, 0.01, 1, 5, 5);
    }

    TEST(cosine_354){
         runTest(3, 3, 0.01, 1, 5, 6);
    }

    TEST(cosine_355){
         runTest(3, 3, 0.01, 1, 5, 7);
    }

    TEST(cosine_356){
         runTest(3, 3, 0.01, 1, 6, 1);
    }

    TEST(cosine_357){
         runTest(3, 3, 0.01, 1, 6, 2);
    }

    TEST(cosine_358){
         runTest(3, 3, 0.01, 1, 6, 3);
    }

    TEST(cosine_359){
         runTest(3, 3, 0.01, 1, 6, 5);
    }

    TEST(cosine_360){
         runTest(3, 3, 0.01, 1, 6, 6);
    }

    TEST(cosine_361){
         runTest(3, 3, 0.01, 1, 6, 7);
    }

    TEST(cosine_362){
         runTest(3, 3, 0.01, 1, 7, 1);
    }

    TEST(cosine_363){
         runTest(3, 3, 0.01, 1, 7, 2);
    }

    TEST(cosine_364){
         runTest(3, 3, 0.01, 1, 7, 3);
    }

    TEST(cosine_365){
         runTest(3, 3, 0.01, 1, 7, 5);
    }

    TEST(cosine_366){
         runTest(3, 3, 0.01, 1, 7, 6);
    }

    TEST(cosine_367){
         runTest(3, 3, 0.01, 1, 7, 7);
    }

    TEST(cosine_368){
         runTest(3, 3, 0.01, 2, 1, 1);
    }

    TEST(cosine_369){
         runTest(3, 3, 0.01, 2, 1, 2);
    }

    TEST(cosine_370){
         runTest(3, 3, 0.01, 2, 1, 3);
    }

    TEST(cosine_371){
         runTest(3, 3, 0.01, 2, 1, 5);
    }

    TEST(cosine_372){
         runTest(3, 3, 0.01, 2, 1, 6);
    }

    TEST(cosine_373){
         runTest(3, 3, 0.01, 2, 1, 7);
    }

    TEST(cosine_374){
         runTest(3, 3, 0.01, 2, 2, 1);
    }

    TEST(cosine_375){
         runTest(3, 3, 0.01, 2, 2, 2);
    }

    TEST(cosine_376){
         runTest(3, 3, 0.01, 2, 2, 3);
    }

    TEST(cosine_377){
         runTest(3, 3, 0.01, 2, 2, 5);
    }

    TEST(cosine_378){
         runTest(3, 3, 0.01, 2, 2, 6);
    }

    TEST(cosine_379){
         runTest(3, 3, 0.01, 2, 2, 7);
    }

    TEST(cosine_380){
         runTest(3, 3, 0.01, 2, 3, 1);
    }

    TEST(cosine_381){
         runTest(3, 3, 0.01, 2, 3, 2);
    }

    TEST(cosine_382){
         runTest(3, 3, 0.01, 2, 3, 3);
    }

    TEST(cosine_383){
         runTest(3, 3, 0.01, 2, 3, 5);
    }

    TEST(cosine_384){
         runTest(3, 3, 0.01, 2, 3, 6);
    }

    TEST(cosine_385){
         runTest(3, 3, 0.01, 2, 3, 7);
    }

    TEST(cosine_386){
         runTest(3, 3, 0.01, 2, 5, 1);
    }

    TEST(cosine_387){
         runTest(3, 3, 0.01, 2, 5, 2);
    }

    TEST(cosine_388){
         runTest(3, 3, 0.01, 2, 5, 3);
    }

    TEST(cosine_389){
         runTest(3, 3, 0.01, 2, 5, 5);
    }

    TEST(cosine_390){
         runTest(3, 3, 0.01, 2, 5, 6);
    }

    TEST(cosine_391){
         runTest(3, 3, 0.01, 2, 5, 7);
    }

    TEST(cosine_392){
         runTest(3, 3, 0.01, 2, 6, 1);
    }

    TEST(cosine_393){
         runTest(3, 3, 0.01, 2, 6, 2);
    }

    TEST(cosine_394){
         runTest(3, 3, 0.01, 2, 6, 3);
    }

    TEST(cosine_395){
         runTest(3, 3, 0.01, 2, 6, 5);
    }

    TEST(cosine_396){
         runTest(3, 3, 0.01, 2, 6, 6);
    }

    TEST(cosine_397){
         runTest(3, 3, 0.01, 2, 6, 7);
    }

    TEST(cosine_398){
         runTest(3, 3, 0.01, 2, 7, 1);
    }

    TEST(cosine_399){
         runTest(3, 3, 0.01, 2, 7, 2);
    }

    TEST(cosine_400){
         runTest(3, 3, 0.01, 2, 7, 3);
    }

    TEST(cosine_401){
         runTest(3, 3, 0.01, 2, 7, 5);
    }

    TEST(cosine_402){
         runTest(3, 3, 0.01, 2, 7, 6);
    }

    TEST(cosine_403){
         runTest(3, 3, 0.01, 2, 7, 7);
    }

    TEST(cosine_404){
         runTest(3, 3, 0.01, 3, 1, 1);
    }

    TEST(cosine_405){
         runTest(3, 3, 0.01, 3, 1, 2);
    }

    TEST(cosine_406){
         runTest(3, 3, 0.01, 3, 1, 3);
    }

    TEST(cosine_407){
         runTest(3, 3, 0.01, 3, 1, 5);
    }

    TEST(cosine_408){
         runTest(3, 3, 0.01, 3, 1, 6);
    }

    TEST(cosine_409){
         runTest(3, 3, 0.01, 3, 1, 7);
    }

    TEST(cosine_410){
         runTest(3, 3, 0.01, 3, 2, 1);
    }

    TEST(cosine_411){
         runTest(3, 3, 0.01, 3, 2, 2);
    }

    TEST(cosine_412){
         runTest(3, 3, 0.01, 3, 2, 3);
    }

    TEST(cosine_413){
         runTest(3, 3, 0.01, 3, 2, 5);
    }

    TEST(cosine_414){
         runTest(3, 3, 0.01, 3, 2, 6);
    }

    TEST(cosine_415){
         runTest(3, 3, 0.01, 3, 2, 7);
    }

    TEST(cosine_416){
         runTest(3, 3, 0.01, 3, 3, 1);
    }

    TEST(cosine_417){
         runTest(3, 3, 0.01, 3, 3, 2);
    }

    TEST(cosine_418){
         runTest(3, 3, 0.01, 3, 3, 3);
    }

    TEST(cosine_419){
         runTest(3, 3, 0.01, 3, 3, 5);
    }

    TEST(cosine_420){
         runTest(3, 3, 0.01, 3, 3, 6);
    }

    TEST(cosine_421){
         runTest(3, 3, 0.01, 3, 3, 7);
    }

    TEST(cosine_422){
         runTest(3, 3, 0.01, 3, 5, 1);
    }

    TEST(cosine_423){
         runTest(3, 3, 0.01, 3, 5, 2);
    }

    TEST(cosine_424){
         runTest(3, 3, 0.01, 3, 5, 3);
    }

    TEST(cosine_425){
         runTest(3, 3, 0.01, 3, 5, 5);
    }

    TEST(cosine_426){
         runTest(3, 3, 0.01, 3, 5, 6);
    }

    TEST(cosine_427){
         runTest(3, 3, 0.01, 3, 5, 7);
    }

    TEST(cosine_428){
         runTest(3, 3, 0.01, 3, 6, 1);
    }

    TEST(cosine_429){
         runTest(3, 3, 0.01, 3, 6, 2);
    }

    TEST(cosine_430){
         runTest(3, 3, 0.01, 3, 6, 3);
    }

    TEST(cosine_431){
         runTest(3, 3, 0.01, 3, 6, 5);
    }

    TEST(cosine_432){
         runTest(3, 3, 0.01, 3, 6, 6);
    }

    TEST(cosine_433){
         runTest(3, 3, 0.01, 3, 6, 7);
    }

    TEST(cosine_434){
         runTest(3, 3, 0.01, 3, 7, 1);
    }

    TEST(cosine_435){
         runTest(3, 3, 0.01, 3, 7, 2);
    }

    TEST(cosine_436){
         runTest(3, 3, 0.01, 3, 7, 3);
    }

    TEST(cosine_437){
         runTest(3, 3, 0.01, 3, 7, 5);
    }

    TEST(cosine_438){
         runTest(3, 3, 0.01, 3, 7, 6);
    }

    TEST(cosine_439){
         runTest(3, 3, 0.01, 3, 7, 7);
    }

    TEST(cosine_440){
         runTest(3, 3, 0.01, 5, 1, 1);
    }

    TEST(cosine_441){
         runTest(3, 3, 0.01, 5, 1, 2);
    }

    TEST(cosine_442){
         runTest(3, 3, 0.01, 5, 1, 3);
    }

    TEST(cosine_443){
         runTest(3, 3, 0.01, 5, 1, 5);
    }

    TEST(cosine_444){
         runTest(3, 3, 0.01, 5, 1, 6);
    }

    TEST(cosine_445){
         runTest(3, 3, 0.01, 5, 1, 7);
    }

    TEST(cosine_446){
         runTest(3, 3, 0.01, 5, 2, 1);
    }

    TEST(cosine_447){
         runTest(3, 3, 0.01, 5, 2, 2);
    }

    TEST(cosine_448){
         runTest(3, 3, 0.01, 5, 2, 3);
    }

    TEST(cosine_449){
         runTest(3, 3, 0.01, 5, 2, 5);
    }

    TEST(cosine_450){
         runTest(3, 3, 0.01, 5, 2, 6);
    }

    TEST(cosine_451){
         runTest(3, 3, 0.01, 5, 2, 7);
    }

    TEST(cosine_452){
         runTest(3, 3, 0.01, 5, 3, 1);
    }

    TEST(cosine_453){
         runTest(3, 3, 0.01, 5, 3, 2);
    }

    TEST(cosine_454){
         runTest(3, 3, 0.01, 5, 3, 3);
    }

    TEST(cosine_455){
         runTest(3, 3, 0.01, 5, 3, 5);
    }

    TEST(cosine_456){
         runTest(3, 3, 0.01, 5, 3, 6);
    }

    TEST(cosine_457){
         runTest(3, 3, 0.01, 5, 3, 7);
    }

    TEST(cosine_458){
         runTest(3, 3, 0.01, 5, 5, 1);
    }

    TEST(cosine_459){
         runTest(3, 3, 0.01, 5, 5, 2);
    }

    TEST(cosine_460){
         runTest(3, 3, 0.01, 5, 5, 3);
    }

    TEST(cosine_461){
         runTest(3, 3, 0.01, 5, 5, 5);
    }

    TEST(cosine_462){
         runTest(3, 3, 0.01, 5, 5, 6);
    }

    TEST(cosine_463){
         runTest(3, 3, 0.01, 5, 5, 7);
    }

    TEST(cosine_464){
         runTest(3, 3, 0.01, 5, 6, 1);
    }

    TEST(cosine_465){
         runTest(3, 3, 0.01, 5, 6, 2);
    }

    TEST(cosine_466){
         runTest(3, 3, 0.01, 5, 6, 3);
    }

    TEST(cosine_467){
         runTest(3, 3, 0.01, 5, 6, 5);
    }

    TEST(cosine_468){
         runTest(3, 3, 0.01, 5, 6, 6);
    }

    TEST(cosine_469){
         runTest(3, 3, 0.01, 5, 6, 7);
    }

    TEST(cosine_470){
         runTest(3, 3, 0.01, 5, 7, 1);
    }

    TEST(cosine_471){
         runTest(3, 3, 0.01, 5, 7, 2);
    }

    TEST(cosine_472){
         runTest(3, 3, 0.01, 5, 7, 3);
    }

    TEST(cosine_473){
         runTest(3, 3, 0.01, 5, 7, 5);
    }

    TEST(cosine_474){
         runTest(3, 3, 0.01, 5, 7, 6);
    }

    TEST(cosine_475){
         runTest(3, 3, 0.01, 5, 7, 7);
    }

    TEST(cosine_476){
         runTest(3, 3, 0.01, 6, 1, 1);
    }

    TEST(cosine_477){
         runTest(3, 3, 0.01, 6, 1, 2);
    }

    TEST(cosine_478){
         runTest(3, 3, 0.01, 6, 1, 3);
    }

    TEST(cosine_479){
         runTest(3, 3, 0.01, 6, 1, 5);
    }

    TEST(cosine_480){
         runTest(3, 3, 0.01, 6, 1, 6);
    }

    TEST(cosine_481){
         runTest(3, 3, 0.01, 6, 1, 7);
    }

    TEST(cosine_482){
         runTest(3, 3, 0.01, 6, 2, 1);
    }

    TEST(cosine_483){
         runTest(3, 3, 0.01, 6, 2, 2);
    }

    TEST(cosine_484){
         runTest(3, 3, 0.01, 6, 2, 3);
    }

    TEST(cosine_485){
         runTest(3, 3, 0.01, 6, 2, 5);
    }

    TEST(cosine_486){
         runTest(3, 3, 0.01, 6, 2, 6);
    }

    TEST(cosine_487){
         runTest(3, 3, 0.01, 6, 2, 7);
    }

    TEST(cosine_488){
         runTest(3, 3, 0.01, 6, 3, 1);
    }

    TEST(cosine_489){
         runTest(3, 3, 0.01, 6, 3, 2);
    }

    TEST(cosine_490){
         runTest(3, 3, 0.01, 6, 3, 3);
    }

    TEST(cosine_491){
         runTest(3, 3, 0.01, 6, 3, 5);
    }

    TEST(cosine_492){
         runTest(3, 3, 0.01, 6, 3, 6);
    }

    TEST(cosine_493){
         runTest(3, 3, 0.01, 6, 3, 7);
    }

    TEST(cosine_494){
         runTest(3, 3, 0.01, 6, 5, 1);
    }

    TEST(cosine_495){
         runTest(3, 3, 0.01, 6, 5, 2);
    }

    TEST(cosine_496){
         runTest(3, 3, 0.01, 6, 5, 3);
    }

    TEST(cosine_497){
         runTest(3, 3, 0.01, 6, 5, 5);
    }

    TEST(cosine_498){
         runTest(3, 3, 0.01, 6, 5, 6);
    }

    TEST(cosine_499){
         runTest(3, 3, 0.01, 6, 5, 7);
    }

    TEST(cosine_500){
         runTest(3, 3, 0.01, 6, 6, 1);
    }

    TEST(cosine_501){
         runTest(3, 3, 0.01, 6, 6, 2);
    }

    TEST(cosine_502){
         runTest(3, 3, 0.01, 6, 6, 3);
    }

    TEST(cosine_503){
         runTest(3, 3, 0.01, 6, 6, 5);
    }

    TEST(cosine_504){
         runTest(3, 3, 0.01, 6, 6, 6);
    }

    TEST(cosine_505){
         runTest(3, 3, 0.01, 6, 6, 7);
    }

    TEST(cosine_506){
         runTest(3, 3, 0.01, 6, 7, 1);
    }

    TEST(cosine_507){
         runTest(3, 3, 0.01, 6, 7, 2);
    }

    TEST(cosine_508){
         runTest(3, 3, 0.01, 6, 7, 3);
    }

    TEST(cosine_509){
         runTest(3, 3, 0.01, 6, 7, 5);
    }

    TEST(cosine_510){
         runTest(3, 3, 0.01, 6, 7, 6);
    }

    TEST(cosine_511){
         runTest(3, 3, 0.01, 6, 7, 7);
    }

    TEST(cosine_512){
         runTest(3, 3, 0.01, 7, 1, 1);
    }

    TEST(cosine_513){
         runTest(3, 3, 0.01, 7, 1, 2);
    }

    TEST(cosine_514){
         runTest(3, 3, 0.01, 7, 1, 3);
    }

    TEST(cosine_515){
         runTest(3, 3, 0.01, 7, 1, 5);
    }

    TEST(cosine_516){
         runTest(3, 3, 0.01, 7, 1, 6);
    }

    TEST(cosine_517){
         runTest(3, 3, 0.01, 7, 1, 7);
    }

    TEST(cosine_518){
         runTest(3, 3, 0.01, 7, 2, 1);
    }

    TEST(cosine_519){
         runTest(3, 3, 0.01, 7, 2, 2);
    }

    TEST(cosine_520){
         runTest(3, 3, 0.01, 7, 2, 3);
    }

    TEST(cosine_521){
         runTest(3, 3, 0.01, 7, 2, 5);
    }

    TEST(cosine_522){
         runTest(3, 3, 0.01, 7, 2, 6);
    }

    TEST(cosine_523){
         runTest(3, 3, 0.01, 7, 2, 7);
    }

    TEST(cosine_524){
         runTest(3, 3, 0.01, 7, 3, 1);
    }

    TEST(cosine_525){
         runTest(3, 3, 0.01, 7, 3, 2);
    }

    TEST(cosine_526){
         runTest(3, 3, 0.01, 7, 3, 3);
    }

    TEST(cosine_527){
         runTest(3, 3, 0.01, 7, 3, 5);
    }

    TEST(cosine_528){
         runTest(3, 3, 0.01, 7, 3, 6);
    }

    TEST(cosine_529){
         runTest(3, 3, 0.01, 7, 3, 7);
    }

    TEST(cosine_530){
         runTest(3, 3, 0.01, 7, 5, 1);
    }

    TEST(cosine_531){
         runTest(3, 3, 0.01, 7, 5, 2);
    }

    TEST(cosine_532){
         runTest(3, 3, 0.01, 7, 5, 3);
    }

    TEST(cosine_533){
         runTest(3, 3, 0.01, 7, 5, 5);
    }

    TEST(cosine_534){
         runTest(3, 3, 0.01, 7, 5, 6);
    }

    TEST(cosine_535){
         runTest(3, 3, 0.01, 7, 5, 7);
    }

    TEST(cosine_536){
         runTest(3, 3, 0.01, 7, 6, 1);
    }

    TEST(cosine_537){
         runTest(3, 3, 0.01, 7, 6, 2);
    }

    TEST(cosine_538){
         runTest(3, 3, 0.01, 7, 6, 3);
    }

    TEST(cosine_539){
         runTest(3, 3, 0.01, 7, 6, 5);
    }

    TEST(cosine_540){
         runTest(3, 3, 0.01, 7, 6, 6);
    }

    TEST(cosine_541){
         runTest(3, 3, 0.01, 7, 6, 7);
    }

    TEST(cosine_542){
         runTest(3, 3, 0.01, 7, 7, 1);
    }

    TEST(cosine_543){
         runTest(3, 3, 0.01, 7, 7, 2);
    }

    TEST(cosine_544){
         runTest(3, 3, 0.01, 7, 7, 3);
    }

    TEST(cosine_545){
         runTest(3, 3, 0.01, 7, 7, 5);
    }

    TEST(cosine_546){
         runTest(3, 3, 0.01, 7, 7, 6);
    }

    TEST(cosine_547){
         runTest(3, 3, 0.01, 7, 7, 7);
    }

    TEST(cosine_548){
         runTest(3, 3, 0.5, 0, 1, 1);
    }

    TEST(cosine_549){
         runTest(3, 3, 0.5, 0, 1, 2);
    }

    TEST(cosine_550){
         runTest(3, 3, 0.5, 0, 1, 3);
    }

    TEST(cosine_551){
         runTest(3, 3, 0.5, 0, 1, 5);
    }

    TEST(cosine_552){
         runTest(3, 3, 0.5, 0, 1, 6);
    }

    TEST(cosine_553){
         runTest(3, 3, 0.5, 0, 1, 7);
    }

    TEST(cosine_554){
         runTest(3, 3, 0.5, 0, 2, 1);
    }

    TEST(cosine_555){
         runTest(3, 3, 0.5, 0, 2, 2);
    }

    TEST(cosine_556){
         runTest(3, 3, 0.5, 0, 2, 3);
    }

    TEST(cosine_557){
         runTest(3, 3, 0.5, 0, 2, 5);
    }

    TEST(cosine_558){
         runTest(3, 3, 0.5, 0, 2, 6);
    }

    TEST(cosine_559){
         runTest(3, 3, 0.5, 0, 2, 7);
    }

    TEST(cosine_560){
         runTest(3, 3, 0.5, 0, 3, 1);
    }

    TEST(cosine_561){
         runTest(3, 3, 0.5, 0, 3, 2);
    }

    TEST(cosine_562){
         runTest(3, 3, 0.5, 0, 3, 3);
    }

    TEST(cosine_563){
         runTest(3, 3, 0.5, 0, 3, 5);
    }

    TEST(cosine_564){
         runTest(3, 3, 0.5, 0, 3, 6);
    }

    TEST(cosine_565){
         runTest(3, 3, 0.5, 0, 3, 7);
    }

    TEST(cosine_566){
         runTest(3, 3, 0.5, 0, 5, 1);
    }

    TEST(cosine_567){
         runTest(3, 3, 0.5, 0, 5, 2);
    }

    TEST(cosine_568){
         runTest(3, 3, 0.5, 0, 5, 3);
    }

    TEST(cosine_569){
         runTest(3, 3, 0.5, 0, 5, 5);
    }

    TEST(cosine_570){
         runTest(3, 3, 0.5, 0, 5, 6);
    }

    TEST(cosine_571){
         runTest(3, 3, 0.5, 0, 5, 7);
    }

    TEST(cosine_572){
         runTest(3, 3, 0.5, 0, 6, 1);
    }

    TEST(cosine_573){
         runTest(3, 3, 0.5, 0, 6, 2);
    }

    TEST(cosine_574){
         runTest(3, 3, 0.5, 0, 6, 3);
    }

    TEST(cosine_575){
         runTest(3, 3, 0.5, 0, 6, 5);
    }

    TEST(cosine_576){
         runTest(3, 3, 0.5, 0, 6, 6);
    }

    TEST(cosine_577){
         runTest(3, 3, 0.5, 0, 6, 7);
    }

    TEST(cosine_578){
         runTest(3, 3, 0.5, 0, 7, 1);
    }

    TEST(cosine_579){
         runTest(3, 3, 0.5, 0, 7, 2);
    }

    TEST(cosine_580){
         runTest(3, 3, 0.5, 0, 7, 3);
    }

    TEST(cosine_581){
         runTest(3, 3, 0.5, 0, 7, 5);
    }

    TEST(cosine_582){
         runTest(3, 3, 0.5, 0, 7, 6);
    }

    TEST(cosine_583){
         runTest(3, 3, 0.5, 0, 7, 7);
    }

    TEST(cosine_584){
         runTest(3, 3, 0.5, 1, 1, 1);
    }

    TEST(cosine_585){
         runTest(3, 3, 0.5, 1, 1, 2);
    }

    TEST(cosine_586){
         runTest(3, 3, 0.5, 1, 1, 3);
    }

    TEST(cosine_587){
         runTest(3, 3, 0.5, 1, 1, 5);
    }

    TEST(cosine_588){
         runTest(3, 3, 0.5, 1, 1, 6);
    }

    TEST(cosine_589){
         runTest(3, 3, 0.5, 1, 1, 7);
    }

    TEST(cosine_590){
         runTest(3, 3, 0.5, 1, 2, 1);
    }

    TEST(cosine_591){
         runTest(3, 3, 0.5, 1, 2, 2);
    }

    TEST(cosine_592){
         runTest(3, 3, 0.5, 1, 2, 3);
    }

    TEST(cosine_593){
         runTest(3, 3, 0.5, 1, 2, 5);
    }

    TEST(cosine_594){
         runTest(3, 3, 0.5, 1, 2, 6);
    }

    TEST(cosine_595){
         runTest(3, 3, 0.5, 1, 2, 7);
    }

    TEST(cosine_596){
         runTest(3, 3, 0.5, 1, 3, 1);
    }

    TEST(cosine_597){
         runTest(3, 3, 0.5, 1, 3, 2);
    }

    TEST(cosine_598){
         runTest(3, 3, 0.5, 1, 3, 3);
    }

    TEST(cosine_599){
         runTest(3, 3, 0.5, 1, 3, 5);
    }

    TEST(cosine_600){
         runTest(3, 3, 0.5, 1, 3, 6);
    }

    TEST(cosine_601){
         runTest(3, 3, 0.5, 1, 3, 7);
    }

    TEST(cosine_602){
         runTest(3, 3, 0.5, 1, 5, 1);
    }

    TEST(cosine_603){
         runTest(3, 3, 0.5, 1, 5, 2);
    }

    TEST(cosine_604){
         runTest(3, 3, 0.5, 1, 5, 3);
    }

    TEST(cosine_605){
         runTest(3, 3, 0.5, 1, 5, 5);
    }

    TEST(cosine_606){
         runTest(3, 3, 0.5, 1, 5, 6);
    }

    TEST(cosine_607){
         runTest(3, 3, 0.5, 1, 5, 7);
    }

    TEST(cosine_608){
         runTest(3, 3, 0.5, 1, 6, 1);
    }

    TEST(cosine_609){
         runTest(3, 3, 0.5, 1, 6, 2);
    }

    TEST(cosine_610){
         runTest(3, 3, 0.5, 1, 6, 3);
    }

    TEST(cosine_611){
         runTest(3, 3, 0.5, 1, 6, 5);
    }

    TEST(cosine_612){
         runTest(3, 3, 0.5, 1, 6, 6);
    }

    TEST(cosine_613){
         runTest(3, 3, 0.5, 1, 6, 7);
    }

    TEST(cosine_614){
         runTest(3, 3, 0.5, 1, 7, 1);
    }

    TEST(cosine_615){
         runTest(3, 3, 0.5, 1, 7, 2);
    }

    TEST(cosine_616){
         runTest(3, 3, 0.5, 1, 7, 3);
    }

    TEST(cosine_617){
         runTest(3, 3, 0.5, 1, 7, 5);
    }

    TEST(cosine_618){
         runTest(3, 3, 0.5, 1, 7, 6);
    }

    TEST(cosine_619){
         runTest(3, 3, 0.5, 1, 7, 7);
    }

    TEST(cosine_620){
         runTest(3, 3, 0.5, 2, 1, 1);
    }

    TEST(cosine_621){
         runTest(3, 3, 0.5, 2, 1, 2);
    }

    TEST(cosine_622){
         runTest(3, 3, 0.5, 2, 1, 3);
    }

    TEST(cosine_623){
         runTest(3, 3, 0.5, 2, 1, 5);
    }

    TEST(cosine_624){
         runTest(3, 3, 0.5, 2, 1, 6);
    }

    TEST(cosine_625){
         runTest(3, 3, 0.5, 2, 1, 7);
    }

    TEST(cosine_626){
         runTest(3, 3, 0.5, 2, 2, 1);
    }

    TEST(cosine_627){
         runTest(3, 3, 0.5, 2, 2, 2);
    }

    TEST(cosine_628){
         runTest(3, 3, 0.5, 2, 2, 3);
    }

    TEST(cosine_629){
         runTest(3, 3, 0.5, 2, 2, 5);
    }

    TEST(cosine_630){
         runTest(3, 3, 0.5, 2, 2, 6);
    }

    TEST(cosine_631){
         runTest(3, 3, 0.5, 2, 2, 7);
    }

    TEST(cosine_632){
         runTest(3, 3, 0.5, 2, 3, 1);
    }

    TEST(cosine_633){
         runTest(3, 3, 0.5, 2, 3, 2);
    }

    TEST(cosine_634){
         runTest(3, 3, 0.5, 2, 3, 3);
    }

    TEST(cosine_635){
         runTest(3, 3, 0.5, 2, 3, 5);
    }

    TEST(cosine_636){
         runTest(3, 3, 0.5, 2, 3, 6);
    }

    TEST(cosine_637){
         runTest(3, 3, 0.5, 2, 3, 7);
    }

    TEST(cosine_638){
         runTest(3, 3, 0.5, 2, 5, 1);
    }

    TEST(cosine_639){
         runTest(3, 3, 0.5, 2, 5, 2);
    }

    TEST(cosine_640){
         runTest(3, 3, 0.5, 2, 5, 3);
    }

    TEST(cosine_641){
         runTest(3, 3, 0.5, 2, 5, 5);
    }

    TEST(cosine_642){
         runTest(3, 3, 0.5, 2, 5, 6);
    }

    TEST(cosine_643){
         runTest(3, 3, 0.5, 2, 5, 7);
    }

    TEST(cosine_644){
         runTest(3, 3, 0.5, 2, 6, 1);
    }

    TEST(cosine_645){
         runTest(3, 3, 0.5, 2, 6, 2);
    }

    TEST(cosine_646){
         runTest(3, 3, 0.5, 2, 6, 3);
    }

    TEST(cosine_647){
         runTest(3, 3, 0.5, 2, 6, 5);
    }

    TEST(cosine_648){
         runTest(3, 3, 0.5, 2, 6, 6);
    }

    TEST(cosine_649){
         runTest(3, 3, 0.5, 2, 6, 7);
    }

    TEST(cosine_650){
         runTest(3, 3, 0.5, 2, 7, 1);
    }

    TEST(cosine_651){
         runTest(3, 3, 0.5, 2, 7, 2);
    }

    TEST(cosine_652){
         runTest(3, 3, 0.5, 2, 7, 3);
    }

    TEST(cosine_653){
         runTest(3, 3, 0.5, 2, 7, 5);
    }

    TEST(cosine_654){
         runTest(3, 3, 0.5, 2, 7, 6);
    }

    TEST(cosine_655){
         runTest(3, 3, 0.5, 2, 7, 7);
    }

    TEST(cosine_656){
         runTest(3, 3, 0.5, 3, 1, 1);
    }

    TEST(cosine_657){
         runTest(3, 3, 0.5, 3, 1, 2);
    }

    TEST(cosine_658){
         runTest(3, 3, 0.5, 3, 1, 3);
    }

    TEST(cosine_659){
         runTest(3, 3, 0.5, 3, 1, 5);
    }

    TEST(cosine_660){
         runTest(3, 3, 0.5, 3, 1, 6);
    }

    TEST(cosine_661){
         runTest(3, 3, 0.5, 3, 1, 7);
    }

    TEST(cosine_662){
         runTest(3, 3, 0.5, 3, 2, 1);
    }

    TEST(cosine_663){
         runTest(3, 3, 0.5, 3, 2, 2);
    }

    TEST(cosine_664){
         runTest(3, 3, 0.5, 3, 2, 3);
    }

    TEST(cosine_665){
         runTest(3, 3, 0.5, 3, 2, 5);
    }

    TEST(cosine_666){
         runTest(3, 3, 0.5, 3, 2, 6);
    }

    TEST(cosine_667){
         runTest(3, 3, 0.5, 3, 2, 7);
    }

    TEST(cosine_668){
         runTest(3, 3, 0.5, 3, 3, 1);
    }

    TEST(cosine_669){
         runTest(3, 3, 0.5, 3, 3, 2);
    }

    TEST(cosine_670){
         runTest(3, 3, 0.5, 3, 3, 3);
    }

    TEST(cosine_671){
         runTest(3, 3, 0.5, 3, 3, 5);
    }

    TEST(cosine_672){
         runTest(3, 3, 0.5, 3, 3, 6);
    }

    TEST(cosine_673){
         runTest(3, 3, 0.5, 3, 3, 7);
    }

    TEST(cosine_674){
         runTest(3, 3, 0.5, 3, 5, 1);
    }

    TEST(cosine_675){
         runTest(3, 3, 0.5, 3, 5, 2);
    }

    TEST(cosine_676){
         runTest(3, 3, 0.5, 3, 5, 3);
    }

    TEST(cosine_677){
         runTest(3, 3, 0.5, 3, 5, 5);
    }

    TEST(cosine_678){
         runTest(3, 3, 0.5, 3, 5, 6);
    }

    TEST(cosine_679){
         runTest(3, 3, 0.5, 3, 5, 7);
    }

    TEST(cosine_680){
         runTest(3, 3, 0.5, 3, 6, 1);
    }

    TEST(cosine_681){
         runTest(3, 3, 0.5, 3, 6, 2);
    }

    TEST(cosine_682){
         runTest(3, 3, 0.5, 3, 6, 3);
    }

    TEST(cosine_683){
         runTest(3, 3, 0.5, 3, 6, 5);
    }

    TEST(cosine_684){
         runTest(3, 3, 0.5, 3, 6, 6);
    }

    TEST(cosine_685){
         runTest(3, 3, 0.5, 3, 6, 7);
    }

    TEST(cosine_686){
         runTest(3, 3, 0.5, 3, 7, 1);
    }

    TEST(cosine_687){
         runTest(3, 3, 0.5, 3, 7, 2);
    }

    TEST(cosine_688){
         runTest(3, 3, 0.5, 3, 7, 3);
    }

    TEST(cosine_689){
         runTest(3, 3, 0.5, 3, 7, 5);
    }

    TEST(cosine_690){
         runTest(3, 3, 0.5, 3, 7, 6);
    }

    TEST(cosine_691){
         runTest(3, 3, 0.5, 3, 7, 7);
    }

    TEST(cosine_692){
         runTest(3, 3, 0.5, 5, 1, 1);
    }

    TEST(cosine_693){
         runTest(3, 3, 0.5, 5, 1, 2);
    }

    TEST(cosine_694){
         runTest(3, 3, 0.5, 5, 1, 3);
    }

    TEST(cosine_695){
         runTest(3, 3, 0.5, 5, 1, 5);
    }

    TEST(cosine_696){
         runTest(3, 3, 0.5, 5, 1, 6);
    }

    TEST(cosine_697){
         runTest(3, 3, 0.5, 5, 1, 7);
    }

    TEST(cosine_698){
         runTest(3, 3, 0.5, 5, 2, 1);
    }

    TEST(cosine_699){
         runTest(3, 3, 0.5, 5, 2, 2);
    }

    TEST(cosine_700){
         runTest(3, 3, 0.5, 5, 2, 3);
    }

    TEST(cosine_701){
         runTest(3, 3, 0.5, 5, 2, 5);
    }

    TEST(cosine_702){
         runTest(3, 3, 0.5, 5, 2, 6);
    }

    TEST(cosine_703){
         runTest(3, 3, 0.5, 5, 2, 7);
    }

    TEST(cosine_704){
         runTest(3, 3, 0.5, 5, 3, 1);
    }

    TEST(cosine_705){
         runTest(3, 3, 0.5, 5, 3, 2);
    }

    TEST(cosine_706){
         runTest(3, 3, 0.5, 5, 3, 3);
    }

    TEST(cosine_707){
         runTest(3, 3, 0.5, 5, 3, 5);
    }

    TEST(cosine_708){
         runTest(3, 3, 0.5, 5, 3, 6);
    }

    TEST(cosine_709){
         runTest(3, 3, 0.5, 5, 3, 7);
    }

    TEST(cosine_710){
         runTest(3, 3, 0.5, 5, 5, 1);
    }

    TEST(cosine_711){
         runTest(3, 3, 0.5, 5, 5, 2);
    }

    TEST(cosine_712){
         runTest(3, 3, 0.5, 5, 5, 3);
    }

    TEST(cosine_713){
         runTest(3, 3, 0.5, 5, 5, 5);
    }

    TEST(cosine_714){
         runTest(3, 3, 0.5, 5, 5, 6);
    }

    TEST(cosine_715){
         runTest(3, 3, 0.5, 5, 5, 7);
    }

    TEST(cosine_716){
         runTest(3, 3, 0.5, 5, 6, 1);
    }

    TEST(cosine_717){
         runTest(3, 3, 0.5, 5, 6, 2);
    }

    TEST(cosine_718){
         runTest(3, 3, 0.5, 5, 6, 3);
    }

    TEST(cosine_719){
         runTest(3, 3, 0.5, 5, 6, 5);
    }

    TEST(cosine_720){
         runTest(3, 3, 0.5, 5, 6, 6);
    }

    TEST(cosine_721){
         runTest(3, 3, 0.5, 5, 6, 7);
    }

    TEST(cosine_722){
         runTest(3, 3, 0.5, 5, 7, 1);
    }

    TEST(cosine_723){
         runTest(3, 3, 0.5, 5, 7, 2);
    }

    TEST(cosine_724){
         runTest(3, 3, 0.5, 5, 7, 3);
    }

    TEST(cosine_725){
         runTest(3, 3, 0.5, 5, 7, 5);
    }

    TEST(cosine_726){
         runTest(3, 3, 0.5, 5, 7, 6);
    }

    TEST(cosine_727){
         runTest(3, 3, 0.5, 5, 7, 7);
    }

    TEST(cosine_728){
         runTest(3, 3, 0.5, 6, 1, 1);
    }

    TEST(cosine_729){
         runTest(3, 3, 0.5, 6, 1, 2);
    }

    TEST(cosine_730){
         runTest(3, 3, 0.5, 6, 1, 3);
    }

    TEST(cosine_731){
         runTest(3, 3, 0.5, 6, 1, 5);
    }

    TEST(cosine_732){
         runTest(3, 3, 0.5, 6, 1, 6);
    }

    TEST(cosine_733){
         runTest(3, 3, 0.5, 6, 1, 7);
    }

    TEST(cosine_734){
         runTest(3, 3, 0.5, 6, 2, 1);
    }

    TEST(cosine_735){
         runTest(3, 3, 0.5, 6, 2, 2);
    }

    TEST(cosine_736){
         runTest(3, 3, 0.5, 6, 2, 3);
    }

    TEST(cosine_737){
         runTest(3, 3, 0.5, 6, 2, 5);
    }

    TEST(cosine_738){
         runTest(3, 3, 0.5, 6, 2, 6);
    }

    TEST(cosine_739){
         runTest(3, 3, 0.5, 6, 2, 7);
    }

    TEST(cosine_740){
         runTest(3, 3, 0.5, 6, 3, 1);
    }

    TEST(cosine_741){
         runTest(3, 3, 0.5, 6, 3, 2);
    }

    TEST(cosine_742){
         runTest(3, 3, 0.5, 6, 3, 3);
    }

    TEST(cosine_743){
         runTest(3, 3, 0.5, 6, 3, 5);
    }

    TEST(cosine_744){
         runTest(3, 3, 0.5, 6, 3, 6);
    }

    TEST(cosine_745){
         runTest(3, 3, 0.5, 6, 3, 7);
    }

    TEST(cosine_746){
         runTest(3, 3, 0.5, 6, 5, 1);
    }

    TEST(cosine_747){
         runTest(3, 3, 0.5, 6, 5, 2);
    }

    TEST(cosine_748){
         runTest(3, 3, 0.5, 6, 5, 3);
    }

    TEST(cosine_749){
         runTest(3, 3, 0.5, 6, 5, 5);
    }

    TEST(cosine_750){
         runTest(3, 3, 0.5, 6, 5, 6);
    }

    TEST(cosine_751){
         runTest(3, 3, 0.5, 6, 5, 7);
    }

    TEST(cosine_752){
         runTest(3, 3, 0.5, 6, 6, 1);
    }

    TEST(cosine_753){
         runTest(3, 3, 0.5, 6, 6, 2);
    }

    TEST(cosine_754){
         runTest(3, 3, 0.5, 6, 6, 3);
    }

    TEST(cosine_755){
         runTest(3, 3, 0.5, 6, 6, 5);
    }

    TEST(cosine_756){
         runTest(3, 3, 0.5, 6, 6, 6);
    }

    TEST(cosine_757){
         runTest(3, 3, 0.5, 6, 6, 7);
    }

    TEST(cosine_758){
         runTest(3, 3, 0.5, 6, 7, 1);
    }

    TEST(cosine_759){
         runTest(3, 3, 0.5, 6, 7, 2);
    }

    TEST(cosine_760){
         runTest(3, 3, 0.5, 6, 7, 3);
    }

    TEST(cosine_761){
         runTest(3, 3, 0.5, 6, 7, 5);
    }

    TEST(cosine_762){
         runTest(3, 3, 0.5, 6, 7, 6);
    }

    TEST(cosine_763){
         runTest(3, 3, 0.5, 6, 7, 7);
    }

    TEST(cosine_764){
         runTest(3, 3, 0.5, 7, 1, 1);
    }

    TEST(cosine_765){
         runTest(3, 3, 0.5, 7, 1, 2);
    }

    TEST(cosine_766){
         runTest(3, 3, 0.5, 7, 1, 3);
    }

    TEST(cosine_767){
         runTest(3, 3, 0.5, 7, 1, 5);
    }

    TEST(cosine_768){
         runTest(3, 3, 0.5, 7, 1, 6);
    }

    TEST(cosine_769){
         runTest(3, 3, 0.5, 7, 1, 7);
    }

    TEST(cosine_770){
         runTest(3, 3, 0.5, 7, 2, 1);
    }

    TEST(cosine_771){
         runTest(3, 3, 0.5, 7, 2, 2);
    }

    TEST(cosine_772){
         runTest(3, 3, 0.5, 7, 2, 3);
    }

    TEST(cosine_773){
         runTest(3, 3, 0.5, 7, 2, 5);
    }

    TEST(cosine_774){
         runTest(3, 3, 0.5, 7, 2, 6);
    }

    TEST(cosine_775){
         runTest(3, 3, 0.5, 7, 2, 7);
    }

    TEST(cosine_776){
         runTest(3, 3, 0.5, 7, 3, 1);
    }

    TEST(cosine_777){
         runTest(3, 3, 0.5, 7, 3, 2);
    }

    TEST(cosine_778){
         runTest(3, 3, 0.5, 7, 3, 3);
    }

    TEST(cosine_779){
         runTest(3, 3, 0.5, 7, 3, 5);
    }

    TEST(cosine_780){
         runTest(3, 3, 0.5, 7, 3, 6);
    }

    TEST(cosine_781){
         runTest(3, 3, 0.5, 7, 3, 7);
    }

    TEST(cosine_782){
         runTest(3, 3, 0.5, 7, 5, 1);
    }

    TEST(cosine_783){
         runTest(3, 3, 0.5, 7, 5, 2);
    }

    TEST(cosine_784){
         runTest(3, 3, 0.5, 7, 5, 3);
    }

    TEST(cosine_785){
         runTest(3, 3, 0.5, 7, 5, 5);
    }

    TEST(cosine_786){
         runTest(3, 3, 0.5, 7, 5, 6);
    }

    TEST(cosine_787){
         runTest(3, 3, 0.5, 7, 5, 7);
    }

    TEST(cosine_788){
         runTest(3, 3, 0.5, 7, 6, 1);
    }

    TEST(cosine_789){
         runTest(3, 3, 0.5, 7, 6, 2);
    }

    TEST(cosine_790){
         runTest(3, 3, 0.5, 7, 6, 3);
    }

    TEST(cosine_791){
         runTest(3, 3, 0.5, 7, 6, 5);
    }

    TEST(cosine_792){
         runTest(3, 3, 0.5, 7, 6, 6);
    }

    TEST(cosine_793){
         runTest(3, 3, 0.5, 7, 6, 7);
    }

    TEST(cosine_794){
         runTest(3, 3, 0.5, 7, 7, 1);
    }

    TEST(cosine_795){
         runTest(3, 3, 0.5, 7, 7, 2);
    }

    TEST(cosine_796){
         runTest(3, 3, 0.5, 7, 7, 3);
    }

    TEST(cosine_797){
         runTest(3, 3, 0.5, 7, 7, 5);
    }

    TEST(cosine_798){
         runTest(3, 3, 0.5, 7, 7, 6);
    }

    TEST(cosine_799){
         runTest(3, 3, 0.5, 7, 7, 7);
    }
}



