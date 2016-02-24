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

void runIntegratorGratingTest(int ns, int nt, double dt, double ds,
                              double C, int wdId, int kxId, int kyId)
{


    Integrator integrator(nt, dt, ns, ds);
    vec k = integrator.spatialFreqVec();
    vec w = integrator.temporalFreqVec();


    double wd = w(wdId);
    vec2 kd = {k(kxId), k(kyId)};

    FullFieldGrating grating(integrator, kd, wd, C);
    grating.computeFourierTransform();
    grating.computeSpatiotemporal();

    cx_cube grating_fft = integrator.backwardFFT(grating.fourierTransform());
    cx_cube diff = grating.spatioTemporal() - grating_fft;

    cube diff_real = abs(real(diff));
    cube diff_imag = abs(imag(diff));

    //Test
    CHECK_CLOSE(diff_real.max(), 0.0, 1e-9);
    CHECK_CLOSE(diff_imag.max(), 0.0, 1e-9);

}


SUITE(integrator){

    TEST(grating_test_0){
         runIntegratorGratingTest(2, 2, 0.01, 0.01, -0.1, 0, 1, 1);
    }

    TEST(grating_test_1){
         runIntegratorGratingTest(2, 2, 0.01, 0.01, 2.2, 0, 1, 1);
    }

    TEST(grating_test_2){
         runIntegratorGratingTest(2, 2, 0.01, 0.01, -0.1, 0, 1, 3);
    }

    TEST(grating_test_3){
         runIntegratorGratingTest(2, 2, 0.01, 0.01, 2.2, 0, 1, 3);
    }

    TEST(grating_test_4){
         runIntegratorGratingTest(2, 2, 0.01, 0.01, -0.1, 0, 3, 1);
    }

    TEST(grating_test_5){
         runIntegratorGratingTest(2, 2, 0.01, 0.01, 2.2, 0, 3, 1);
    }

    TEST(grating_test_6){
         runIntegratorGratingTest(2, 2, 0.01, 0.01, -0.1, 0, 3, 3);
    }

    TEST(grating_test_7){
         runIntegratorGratingTest(2, 2, 0.01, 0.01, 2.2, 0, 3, 3);
    }

    TEST(grating_test_8){
         runIntegratorGratingTest(2, 2, 0.01, 0.01, -0.1, 1, 1, 1);
    }

    TEST(grating_test_9){
         runIntegratorGratingTest(2, 2, 0.01, 0.01, 2.2, 1, 1, 1);
    }

    TEST(grating_test_10){
         runIntegratorGratingTest(2, 2, 0.01, 0.01, -0.1, 1, 1, 3);
    }

    TEST(grating_test_11){
         runIntegratorGratingTest(2, 2, 0.01, 0.01, 2.2, 1, 1, 3);
    }

    TEST(grating_test_12){
         runIntegratorGratingTest(2, 2, 0.01, 0.01, -0.1, 1, 3, 1);
    }

    TEST(grating_test_13){
         runIntegratorGratingTest(2, 2, 0.01, 0.01, 2.2, 1, 3, 1);
    }

    TEST(grating_test_14){
         runIntegratorGratingTest(2, 2, 0.01, 0.01, -0.1, 1, 3, 3);
    }

    TEST(grating_test_15){
         runIntegratorGratingTest(2, 2, 0.01, 0.01, 2.2, 1, 3, 3);
    }

    TEST(grating_test_16){
         runIntegratorGratingTest(2, 2, 0.01, 0.01, -0.1, 3, 1, 1);
    }

    TEST(grating_test_17){
         runIntegratorGratingTest(2, 2, 0.01, 0.01, 2.2, 3, 1, 1);
    }

    TEST(grating_test_18){
         runIntegratorGratingTest(2, 2, 0.01, 0.01, -0.1, 3, 1, 3);
    }

    TEST(grating_test_19){
         runIntegratorGratingTest(2, 2, 0.01, 0.01, 2.2, 3, 1, 3);
    }

    TEST(grating_test_20){
         runIntegratorGratingTest(2, 2, 0.01, 0.01, -0.1, 3, 3, 1);
    }

    TEST(grating_test_21){
         runIntegratorGratingTest(2, 2, 0.01, 0.01, 2.2, 3, 3, 1);
    }

    TEST(grating_test_22){
         runIntegratorGratingTest(2, 2, 0.01, 0.01, -0.1, 3, 3, 3);
    }

    TEST(grating_test_23){
         runIntegratorGratingTest(2, 2, 0.01, 0.01, 2.2, 3, 3, 3);
    }

    TEST(grating_test_24){
         runIntegratorGratingTest(2, 2, 0.5, 0.01, -0.1, 0, 1, 1);
    }

    TEST(grating_test_25){
         runIntegratorGratingTest(2, 2, 0.5, 0.01, 2.2, 0, 1, 1);
    }

    TEST(grating_test_26){
         runIntegratorGratingTest(2, 2, 0.5, 0.01, -0.1, 0, 1, 3);
    }

    TEST(grating_test_27){
         runIntegratorGratingTest(2, 2, 0.5, 0.01, 2.2, 0, 1, 3);
    }

    TEST(grating_test_28){
         runIntegratorGratingTest(2, 2, 0.5, 0.01, -0.1, 0, 3, 1);
    }

    TEST(grating_test_29){
         runIntegratorGratingTest(2, 2, 0.5, 0.01, 2.2, 0, 3, 1);
    }

    TEST(grating_test_30){
         runIntegratorGratingTest(2, 2, 0.5, 0.01, -0.1, 0, 3, 3);
    }

    TEST(grating_test_31){
         runIntegratorGratingTest(2, 2, 0.5, 0.01, 2.2, 0, 3, 3);
    }

    TEST(grating_test_32){
         runIntegratorGratingTest(2, 2, 0.5, 0.01, -0.1, 1, 1, 1);
    }

    TEST(grating_test_33){
         runIntegratorGratingTest(2, 2, 0.5, 0.01, 2.2, 1, 1, 1);
    }

    TEST(grating_test_34){
         runIntegratorGratingTest(2, 2, 0.5, 0.01, -0.1, 1, 1, 3);
    }

    TEST(grating_test_35){
         runIntegratorGratingTest(2, 2, 0.5, 0.01, 2.2, 1, 1, 3);
    }

    TEST(grating_test_36){
         runIntegratorGratingTest(2, 2, 0.5, 0.01, -0.1, 1, 3, 1);
    }

    TEST(grating_test_37){
         runIntegratorGratingTest(2, 2, 0.5, 0.01, 2.2, 1, 3, 1);
    }

    TEST(grating_test_38){
         runIntegratorGratingTest(2, 2, 0.5, 0.01, -0.1, 1, 3, 3);
    }

    TEST(grating_test_39){
         runIntegratorGratingTest(2, 2, 0.5, 0.01, 2.2, 1, 3, 3);
    }

    TEST(grating_test_40){
         runIntegratorGratingTest(2, 2, 0.5, 0.01, -0.1, 3, 1, 1);
    }

    TEST(grating_test_41){
         runIntegratorGratingTest(2, 2, 0.5, 0.01, 2.2, 3, 1, 1);
    }

    TEST(grating_test_42){
         runIntegratorGratingTest(2, 2, 0.5, 0.01, -0.1, 3, 1, 3);
    }

    TEST(grating_test_43){
         runIntegratorGratingTest(2, 2, 0.5, 0.01, 2.2, 3, 1, 3);
    }

    TEST(grating_test_44){
         runIntegratorGratingTest(2, 2, 0.5, 0.01, -0.1, 3, 3, 1);
    }

    TEST(grating_test_45){
         runIntegratorGratingTest(2, 2, 0.5, 0.01, 2.2, 3, 3, 1);
    }

    TEST(grating_test_46){
         runIntegratorGratingTest(2, 2, 0.5, 0.01, -0.1, 3, 3, 3);
    }

    TEST(grating_test_47){
         runIntegratorGratingTest(2, 2, 0.5, 0.01, 2.2, 3, 3, 3);
    }

    TEST(grating_test_48){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, -0.1, 0, 1, 1);
    }

    TEST(grating_test_49){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, 2.2, 0, 1, 1);
    }

    TEST(grating_test_50){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, -0.1, 0, 1, 3);
    }

    TEST(grating_test_51){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, 2.2, 0, 1, 3);
    }

    TEST(grating_test_52){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, -0.1, 0, 3, 1);
    }

    TEST(grating_test_53){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, 2.2, 0, 3, 1);
    }

    TEST(grating_test_54){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, -0.1, 0, 3, 3);
    }

    TEST(grating_test_55){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, 2.2, 0, 3, 3);
    }

    TEST(grating_test_56){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, -0.1, 1, 1, 1);
    }

    TEST(grating_test_57){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, 2.2, 1, 1, 1);
    }

    TEST(grating_test_58){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, -0.1, 1, 1, 3);
    }

    TEST(grating_test_59){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, 2.2, 1, 1, 3);
    }

    TEST(grating_test_60){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, -0.1, 1, 3, 1);
    }

    TEST(grating_test_61){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, 2.2, 1, 3, 1);
    }

    TEST(grating_test_62){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, -0.1, 1, 3, 3);
    }

    TEST(grating_test_63){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, 2.2, 1, 3, 3);
    }

    TEST(grating_test_64){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, -0.1, 2, 1, 1);
    }

    TEST(grating_test_65){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, 2.2, 2, 1, 1);
    }

    TEST(grating_test_66){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, -0.1, 2, 1, 3);
    }

    TEST(grating_test_67){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, 2.2, 2, 1, 3);
    }

    TEST(grating_test_68){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, -0.1, 2, 3, 1);
    }

    TEST(grating_test_69){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, 2.2, 2, 3, 1);
    }

    TEST(grating_test_70){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, -0.1, 2, 3, 3);
    }

    TEST(grating_test_71){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, 2.2, 2, 3, 3);
    }

    TEST(grating_test_72){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, -0.1, 3, 1, 1);
    }

    TEST(grating_test_73){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, 2.2, 3, 1, 1);
    }

    TEST(grating_test_74){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, -0.1, 3, 1, 3);
    }

    TEST(grating_test_75){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, 2.2, 3, 1, 3);
    }

    TEST(grating_test_76){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, -0.1, 3, 3, 1);
    }

    TEST(grating_test_77){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, 2.2, 3, 3, 1);
    }

    TEST(grating_test_78){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, -0.1, 3, 3, 3);
    }

    TEST(grating_test_79){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, 2.2, 3, 3, 3);
    }

    TEST(grating_test_80){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, -0.1, 5, 1, 1);
    }

    TEST(grating_test_81){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, 2.2, 5, 1, 1);
    }

    TEST(grating_test_82){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, -0.1, 5, 1, 3);
    }

    TEST(grating_test_83){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, 2.2, 5, 1, 3);
    }

    TEST(grating_test_84){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, -0.1, 5, 3, 1);
    }

    TEST(grating_test_85){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, 2.2, 5, 3, 1);
    }

    TEST(grating_test_86){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, -0.1, 5, 3, 3);
    }

    TEST(grating_test_87){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, 2.2, 5, 3, 3);
    }

    TEST(grating_test_88){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, -0.1, 6, 1, 1);
    }

    TEST(grating_test_89){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, 2.2, 6, 1, 1);
    }

    TEST(grating_test_90){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, -0.1, 6, 1, 3);
    }

    TEST(grating_test_91){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, 2.2, 6, 1, 3);
    }

    TEST(grating_test_92){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, -0.1, 6, 3, 1);
    }

    TEST(grating_test_93){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, 2.2, 6, 3, 1);
    }

    TEST(grating_test_94){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, -0.1, 6, 3, 3);
    }

    TEST(grating_test_95){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, 2.2, 6, 3, 3);
    }

    TEST(grating_test_96){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, -0.1, 7, 1, 1);
    }

    TEST(grating_test_97){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, 2.2, 7, 1, 1);
    }

    TEST(grating_test_98){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, -0.1, 7, 1, 3);
    }

    TEST(grating_test_99){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, 2.2, 7, 1, 3);
    }

    TEST(grating_test_100){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, -0.1, 7, 3, 1);
    }

    TEST(grating_test_101){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, 2.2, 7, 3, 1);
    }

    TEST(grating_test_102){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, -0.1, 7, 3, 3);
    }

    TEST(grating_test_103){
         runIntegratorGratingTest(2, 3, 0.01, 0.01, 2.2, 7, 3, 3);
    }

    TEST(grating_test_104){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, -0.1, 0, 1, 1);
    }

    TEST(grating_test_105){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, 2.2, 0, 1, 1);
    }

    TEST(grating_test_106){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, -0.1, 0, 1, 3);
    }

    TEST(grating_test_107){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, 2.2, 0, 1, 3);
    }

    TEST(grating_test_108){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, -0.1, 0, 3, 1);
    }

    TEST(grating_test_109){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, 2.2, 0, 3, 1);
    }

    TEST(grating_test_110){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, -0.1, 0, 3, 3);
    }

    TEST(grating_test_111){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, 2.2, 0, 3, 3);
    }

    TEST(grating_test_112){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, -0.1, 1, 1, 1);
    }

    TEST(grating_test_113){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, 2.2, 1, 1, 1);
    }

    TEST(grating_test_114){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, -0.1, 1, 1, 3);
    }

    TEST(grating_test_115){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, 2.2, 1, 1, 3);
    }

    TEST(grating_test_116){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, -0.1, 1, 3, 1);
    }

    TEST(grating_test_117){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, 2.2, 1, 3, 1);
    }

    TEST(grating_test_118){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, -0.1, 1, 3, 3);
    }

    TEST(grating_test_119){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, 2.2, 1, 3, 3);
    }

    TEST(grating_test_120){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, -0.1, 2, 1, 1);
    }

    TEST(grating_test_121){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, 2.2, 2, 1, 1);
    }

    TEST(grating_test_122){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, -0.1, 2, 1, 3);
    }

    TEST(grating_test_123){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, 2.2, 2, 1, 3);
    }

    TEST(grating_test_124){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, -0.1, 2, 3, 1);
    }

    TEST(grating_test_125){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, 2.2, 2, 3, 1);
    }

    TEST(grating_test_126){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, -0.1, 2, 3, 3);
    }

    TEST(grating_test_127){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, 2.2, 2, 3, 3);
    }

    TEST(grating_test_128){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, -0.1, 3, 1, 1);
    }

    TEST(grating_test_129){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, 2.2, 3, 1, 1);
    }

    TEST(grating_test_130){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, -0.1, 3, 1, 3);
    }

    TEST(grating_test_131){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, 2.2, 3, 1, 3);
    }

    TEST(grating_test_132){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, -0.1, 3, 3, 1);
    }

    TEST(grating_test_133){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, 2.2, 3, 3, 1);
    }

    TEST(grating_test_134){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, -0.1, 3, 3, 3);
    }

    TEST(grating_test_135){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, 2.2, 3, 3, 3);
    }

    TEST(grating_test_136){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, -0.1, 5, 1, 1);
    }

    TEST(grating_test_137){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, 2.2, 5, 1, 1);
    }

    TEST(grating_test_138){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, -0.1, 5, 1, 3);
    }

    TEST(grating_test_139){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, 2.2, 5, 1, 3);
    }

    TEST(grating_test_140){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, -0.1, 5, 3, 1);
    }

    TEST(grating_test_141){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, 2.2, 5, 3, 1);
    }

    TEST(grating_test_142){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, -0.1, 5, 3, 3);
    }

    TEST(grating_test_143){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, 2.2, 5, 3, 3);
    }

    TEST(grating_test_144){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, -0.1, 6, 1, 1);
    }

    TEST(grating_test_145){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, 2.2, 6, 1, 1);
    }

    TEST(grating_test_146){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, -0.1, 6, 1, 3);
    }

    TEST(grating_test_147){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, 2.2, 6, 1, 3);
    }

    TEST(grating_test_148){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, -0.1, 6, 3, 1);
    }

    TEST(grating_test_149){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, 2.2, 6, 3, 1);
    }

    TEST(grating_test_150){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, -0.1, 6, 3, 3);
    }

    TEST(grating_test_151){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, 2.2, 6, 3, 3);
    }

    TEST(grating_test_152){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, -0.1, 7, 1, 1);
    }

    TEST(grating_test_153){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, 2.2, 7, 1, 1);
    }

    TEST(grating_test_154){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, -0.1, 7, 1, 3);
    }

    TEST(grating_test_155){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, 2.2, 7, 1, 3);
    }

    TEST(grating_test_156){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, -0.1, 7, 3, 1);
    }

    TEST(grating_test_157){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, 2.2, 7, 3, 1);
    }

    TEST(grating_test_158){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, -0.1, 7, 3, 3);
    }

    TEST(grating_test_159){
         runIntegratorGratingTest(2, 3, 0.5, 0.01, 2.2, 7, 3, 3);
    }

    TEST(grating_test_160){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 1, 1);
    }

    TEST(grating_test_161){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 1, 1);
    }

    TEST(grating_test_162){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 1, 2);
    }

    TEST(grating_test_163){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 1, 2);
    }

    TEST(grating_test_164){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 1, 3);
    }

    TEST(grating_test_165){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 1, 3);
    }

    TEST(grating_test_166){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 1, 5);
    }

    TEST(grating_test_167){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 1, 5);
    }

    TEST(grating_test_168){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 1, 6);
    }

    TEST(grating_test_169){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 1, 6);
    }

    TEST(grating_test_170){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 1, 7);
    }

    TEST(grating_test_171){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 1, 7);
    }

    TEST(grating_test_172){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 2, 1);
    }

    TEST(grating_test_173){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 2, 1);
    }

    TEST(grating_test_174){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 2, 2);
    }

    TEST(grating_test_175){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 2, 2);
    }

    TEST(grating_test_176){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 2, 3);
    }

    TEST(grating_test_177){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 2, 3);
    }

    TEST(grating_test_178){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 2, 5);
    }

    TEST(grating_test_179){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 2, 5);
    }

    TEST(grating_test_180){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 2, 6);
    }

    TEST(grating_test_181){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 2, 6);
    }

    TEST(grating_test_182){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 2, 7);
    }

    TEST(grating_test_183){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 2, 7);
    }

    TEST(grating_test_184){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 3, 1);
    }

    TEST(grating_test_185){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 3, 1);
    }

    TEST(grating_test_186){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 3, 2);
    }

    TEST(grating_test_187){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 3, 2);
    }

    TEST(grating_test_188){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 3, 3);
    }

    TEST(grating_test_189){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 3, 3);
    }

    TEST(grating_test_190){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 3, 5);
    }

    TEST(grating_test_191){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 3, 5);
    }

    TEST(grating_test_192){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 3, 6);
    }

    TEST(grating_test_193){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 3, 6);
    }

    TEST(grating_test_194){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 3, 7);
    }

    TEST(grating_test_195){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 3, 7);
    }

    TEST(grating_test_196){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 5, 1);
    }

    TEST(grating_test_197){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 5, 1);
    }

    TEST(grating_test_198){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 5, 2);
    }

    TEST(grating_test_199){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 5, 2);
    }

    TEST(grating_test_200){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 5, 3);
    }

    TEST(grating_test_201){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 5, 3);
    }

    TEST(grating_test_202){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 5, 5);
    }

    TEST(grating_test_203){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 5, 5);
    }

    TEST(grating_test_204){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 5, 6);
    }

    TEST(grating_test_205){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 5, 6);
    }

    TEST(grating_test_206){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 5, 7);
    }

    TEST(grating_test_207){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 5, 7);
    }

    TEST(grating_test_208){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 6, 1);
    }

    TEST(grating_test_209){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 6, 1);
    }

    TEST(grating_test_210){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 6, 2);
    }

    TEST(grating_test_211){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 6, 2);
    }

    TEST(grating_test_212){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 6, 3);
    }

    TEST(grating_test_213){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 6, 3);
    }

    TEST(grating_test_214){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 6, 5);
    }

    TEST(grating_test_215){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 6, 5);
    }

    TEST(grating_test_216){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 6, 6);
    }

    TEST(grating_test_217){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 6, 6);
    }

    TEST(grating_test_218){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 6, 7);
    }

    TEST(grating_test_219){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 6, 7);
    }

    TEST(grating_test_220){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 7, 1);
    }

    TEST(grating_test_221){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 7, 1);
    }

    TEST(grating_test_222){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 7, 2);
    }

    TEST(grating_test_223){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 7, 2);
    }

    TEST(grating_test_224){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 7, 3);
    }

    TEST(grating_test_225){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 7, 3);
    }

    TEST(grating_test_226){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 7, 5);
    }

    TEST(grating_test_227){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 7, 5);
    }

    TEST(grating_test_228){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 7, 6);
    }

    TEST(grating_test_229){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 7, 6);
    }

    TEST(grating_test_230){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 7, 7);
    }

    TEST(grating_test_231){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 7, 7);
    }

    TEST(grating_test_232){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 1, 1);
    }

    TEST(grating_test_233){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 1, 1);
    }

    TEST(grating_test_234){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 1, 2);
    }

    TEST(grating_test_235){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 1, 2);
    }

    TEST(grating_test_236){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 1, 3);
    }

    TEST(grating_test_237){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 1, 3);
    }

    TEST(grating_test_238){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 1, 5);
    }

    TEST(grating_test_239){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 1, 5);
    }

    TEST(grating_test_240){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 1, 6);
    }

    TEST(grating_test_241){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 1, 6);
    }

    TEST(grating_test_242){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 1, 7);
    }

    TEST(grating_test_243){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 1, 7);
    }

    TEST(grating_test_244){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 2, 1);
    }

    TEST(grating_test_245){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 2, 1);
    }

    TEST(grating_test_246){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 2, 2);
    }

    TEST(grating_test_247){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 2, 2);
    }

    TEST(grating_test_248){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 2, 3);
    }

    TEST(grating_test_249){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 2, 3);
    }

    TEST(grating_test_250){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 2, 5);
    }

    TEST(grating_test_251){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 2, 5);
    }

    TEST(grating_test_252){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 2, 6);
    }

    TEST(grating_test_253){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 2, 6);
    }

    TEST(grating_test_254){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 2, 7);
    }

    TEST(grating_test_255){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 2, 7);
    }

    TEST(grating_test_256){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 3, 1);
    }

    TEST(grating_test_257){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 3, 1);
    }

    TEST(grating_test_258){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 3, 2);
    }

    TEST(grating_test_259){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 3, 2);
    }

    TEST(grating_test_260){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 3, 3);
    }

    TEST(grating_test_261){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 3, 3);
    }

    TEST(grating_test_262){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 3, 5);
    }

    TEST(grating_test_263){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 3, 5);
    }

    TEST(grating_test_264){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 3, 6);
    }

    TEST(grating_test_265){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 3, 6);
    }

    TEST(grating_test_266){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 3, 7);
    }

    TEST(grating_test_267){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 3, 7);
    }

    TEST(grating_test_268){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 5, 1);
    }

    TEST(grating_test_269){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 5, 1);
    }

    TEST(grating_test_270){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 5, 2);
    }

    TEST(grating_test_271){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 5, 2);
    }

    TEST(grating_test_272){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 5, 3);
    }

    TEST(grating_test_273){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 5, 3);
    }

    TEST(grating_test_274){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 5, 5);
    }

    TEST(grating_test_275){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 5, 5);
    }

    TEST(grating_test_276){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 5, 6);
    }

    TEST(grating_test_277){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 5, 6);
    }

    TEST(grating_test_278){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 5, 7);
    }

    TEST(grating_test_279){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 5, 7);
    }

    TEST(grating_test_280){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 6, 1);
    }

    TEST(grating_test_281){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 6, 1);
    }

    TEST(grating_test_282){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 6, 2);
    }

    TEST(grating_test_283){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 6, 2);
    }

    TEST(grating_test_284){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 6, 3);
    }

    TEST(grating_test_285){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 6, 3);
    }

    TEST(grating_test_286){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 6, 5);
    }

    TEST(grating_test_287){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 6, 5);
    }

    TEST(grating_test_288){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 6, 6);
    }

    TEST(grating_test_289){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 6, 6);
    }

    TEST(grating_test_290){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 6, 7);
    }

    TEST(grating_test_291){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 6, 7);
    }

    TEST(grating_test_292){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 7, 1);
    }

    TEST(grating_test_293){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 7, 1);
    }

    TEST(grating_test_294){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 7, 2);
    }

    TEST(grating_test_295){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 7, 2);
    }

    TEST(grating_test_296){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 7, 3);
    }

    TEST(grating_test_297){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 7, 3);
    }

    TEST(grating_test_298){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 7, 5);
    }

    TEST(grating_test_299){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 7, 5);
    }

    TEST(grating_test_300){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 7, 6);
    }

    TEST(grating_test_301){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 7, 6);
    }

    TEST(grating_test_302){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 7, 7);
    }

    TEST(grating_test_303){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 7, 7);
    }

    TEST(grating_test_304){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 1, 1);
    }

    TEST(grating_test_305){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 1, 1);
    }

    TEST(grating_test_306){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 1, 2);
    }

    TEST(grating_test_307){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 1, 2);
    }

    TEST(grating_test_308){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 1, 3);
    }

    TEST(grating_test_309){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 1, 3);
    }

    TEST(grating_test_310){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 1, 5);
    }

    TEST(grating_test_311){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 1, 5);
    }

    TEST(grating_test_312){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 1, 6);
    }

    TEST(grating_test_313){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 1, 6);
    }

    TEST(grating_test_314){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 1, 7);
    }

    TEST(grating_test_315){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 1, 7);
    }

    TEST(grating_test_316){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 2, 1);
    }

    TEST(grating_test_317){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 2, 1);
    }

    TEST(grating_test_318){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 2, 2);
    }

    TEST(grating_test_319){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 2, 2);
    }

    TEST(grating_test_320){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 2, 3);
    }

    TEST(grating_test_321){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 2, 3);
    }

    TEST(grating_test_322){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 2, 5);
    }

    TEST(grating_test_323){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 2, 5);
    }

    TEST(grating_test_324){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 2, 6);
    }

    TEST(grating_test_325){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 2, 6);
    }

    TEST(grating_test_326){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 2, 7);
    }

    TEST(grating_test_327){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 2, 7);
    }

    TEST(grating_test_328){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 3, 1);
    }

    TEST(grating_test_329){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 3, 1);
    }

    TEST(grating_test_330){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 3, 2);
    }

    TEST(grating_test_331){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 3, 2);
    }

    TEST(grating_test_332){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 3, 3);
    }

    TEST(grating_test_333){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 3, 3);
    }

    TEST(grating_test_334){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 3, 5);
    }

    TEST(grating_test_335){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 3, 5);
    }

    TEST(grating_test_336){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 3, 6);
    }

    TEST(grating_test_337){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 3, 6);
    }

    TEST(grating_test_338){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 3, 7);
    }

    TEST(grating_test_339){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 3, 7);
    }

    TEST(grating_test_340){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 5, 1);
    }

    TEST(grating_test_341){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 5, 1);
    }

    TEST(grating_test_342){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 5, 2);
    }

    TEST(grating_test_343){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 5, 2);
    }

    TEST(grating_test_344){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 5, 3);
    }

    TEST(grating_test_345){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 5, 3);
    }

    TEST(grating_test_346){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 5, 5);
    }

    TEST(grating_test_347){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 5, 5);
    }

    TEST(grating_test_348){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 5, 6);
    }

    TEST(grating_test_349){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 5, 6);
    }

    TEST(grating_test_350){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 5, 7);
    }

    TEST(grating_test_351){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 5, 7);
    }

    TEST(grating_test_352){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 6, 1);
    }

    TEST(grating_test_353){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 6, 1);
    }

    TEST(grating_test_354){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 6, 2);
    }

    TEST(grating_test_355){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 6, 2);
    }

    TEST(grating_test_356){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 6, 3);
    }

    TEST(grating_test_357){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 6, 3);
    }

    TEST(grating_test_358){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 6, 5);
    }

    TEST(grating_test_359){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 6, 5);
    }

    TEST(grating_test_360){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 6, 6);
    }

    TEST(grating_test_361){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 6, 6);
    }

    TEST(grating_test_362){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 6, 7);
    }

    TEST(grating_test_363){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 6, 7);
    }

    TEST(grating_test_364){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 7, 1);
    }

    TEST(grating_test_365){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 7, 1);
    }

    TEST(grating_test_366){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 7, 2);
    }

    TEST(grating_test_367){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 7, 2);
    }

    TEST(grating_test_368){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 7, 3);
    }

    TEST(grating_test_369){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 7, 3);
    }

    TEST(grating_test_370){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 7, 5);
    }

    TEST(grating_test_371){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 7, 5);
    }

    TEST(grating_test_372){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 7, 6);
    }

    TEST(grating_test_373){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 7, 6);
    }

    TEST(grating_test_374){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 7, 7);
    }

    TEST(grating_test_375){
         runIntegratorGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 7, 7);
    }

    TEST(grating_test_376){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 1, 1);
    }

    TEST(grating_test_377){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 1, 1);
    }

    TEST(grating_test_378){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 1, 2);
    }

    TEST(grating_test_379){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 1, 2);
    }

    TEST(grating_test_380){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 1, 3);
    }

    TEST(grating_test_381){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 1, 3);
    }

    TEST(grating_test_382){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 1, 5);
    }

    TEST(grating_test_383){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 1, 5);
    }

    TEST(grating_test_384){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 1, 6);
    }

    TEST(grating_test_385){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 1, 6);
    }

    TEST(grating_test_386){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 1, 7);
    }

    TEST(grating_test_387){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 1, 7);
    }

    TEST(grating_test_388){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 2, 1);
    }

    TEST(grating_test_389){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 2, 1);
    }

    TEST(grating_test_390){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 2, 2);
    }

    TEST(grating_test_391){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 2, 2);
    }

    TEST(grating_test_392){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 2, 3);
    }

    TEST(grating_test_393){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 2, 3);
    }

    TEST(grating_test_394){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 2, 5);
    }

    TEST(grating_test_395){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 2, 5);
    }

    TEST(grating_test_396){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 2, 6);
    }

    TEST(grating_test_397){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 2, 6);
    }

    TEST(grating_test_398){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 2, 7);
    }

    TEST(grating_test_399){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 2, 7);
    }

    TEST(grating_test_400){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 3, 1);
    }

    TEST(grating_test_401){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 3, 1);
    }

    TEST(grating_test_402){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 3, 2);
    }

    TEST(grating_test_403){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 3, 2);
    }

    TEST(grating_test_404){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 3, 3);
    }

    TEST(grating_test_405){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 3, 3);
    }

    TEST(grating_test_406){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 3, 5);
    }

    TEST(grating_test_407){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 3, 5);
    }

    TEST(grating_test_408){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 3, 6);
    }

    TEST(grating_test_409){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 3, 6);
    }

    TEST(grating_test_410){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 3, 7);
    }

    TEST(grating_test_411){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 3, 7);
    }

    TEST(grating_test_412){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 5, 1);
    }

    TEST(grating_test_413){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 5, 1);
    }

    TEST(grating_test_414){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 5, 2);
    }

    TEST(grating_test_415){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 5, 2);
    }

    TEST(grating_test_416){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 5, 3);
    }

    TEST(grating_test_417){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 5, 3);
    }

    TEST(grating_test_418){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 5, 5);
    }

    TEST(grating_test_419){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 5, 5);
    }

    TEST(grating_test_420){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 5, 6);
    }

    TEST(grating_test_421){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 5, 6);
    }

    TEST(grating_test_422){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 5, 7);
    }

    TEST(grating_test_423){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 5, 7);
    }

    TEST(grating_test_424){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 6, 1);
    }

    TEST(grating_test_425){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 6, 1);
    }

    TEST(grating_test_426){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 6, 2);
    }

    TEST(grating_test_427){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 6, 2);
    }

    TEST(grating_test_428){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 6, 3);
    }

    TEST(grating_test_429){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 6, 3);
    }

    TEST(grating_test_430){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 6, 5);
    }

    TEST(grating_test_431){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 6, 5);
    }

    TEST(grating_test_432){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 6, 6);
    }

    TEST(grating_test_433){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 6, 6);
    }

    TEST(grating_test_434){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 6, 7);
    }

    TEST(grating_test_435){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 6, 7);
    }

    TEST(grating_test_436){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 7, 1);
    }

    TEST(grating_test_437){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 7, 1);
    }

    TEST(grating_test_438){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 7, 2);
    }

    TEST(grating_test_439){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 7, 2);
    }

    TEST(grating_test_440){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 7, 3);
    }

    TEST(grating_test_441){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 7, 3);
    }

    TEST(grating_test_442){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 7, 5);
    }

    TEST(grating_test_443){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 7, 5);
    }

    TEST(grating_test_444){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 7, 6);
    }

    TEST(grating_test_445){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 7, 6);
    }

    TEST(grating_test_446){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 7, 7);
    }

    TEST(grating_test_447){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 7, 7);
    }

    TEST(grating_test_448){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 1, 1);
    }

    TEST(grating_test_449){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 1, 1);
    }

    TEST(grating_test_450){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 1, 2);
    }

    TEST(grating_test_451){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 1, 2);
    }

    TEST(grating_test_452){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 1, 3);
    }

    TEST(grating_test_453){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 1, 3);
    }

    TEST(grating_test_454){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 1, 5);
    }

    TEST(grating_test_455){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 1, 5);
    }

    TEST(grating_test_456){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 1, 6);
    }

    TEST(grating_test_457){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 1, 6);
    }

    TEST(grating_test_458){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 1, 7);
    }

    TEST(grating_test_459){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 1, 7);
    }

    TEST(grating_test_460){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 2, 1);
    }

    TEST(grating_test_461){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 2, 1);
    }

    TEST(grating_test_462){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 2, 2);
    }

    TEST(grating_test_463){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 2, 2);
    }

    TEST(grating_test_464){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 2, 3);
    }

    TEST(grating_test_465){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 2, 3);
    }

    TEST(grating_test_466){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 2, 5);
    }

    TEST(grating_test_467){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 2, 5);
    }

    TEST(grating_test_468){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 2, 6);
    }

    TEST(grating_test_469){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 2, 6);
    }

    TEST(grating_test_470){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 2, 7);
    }

    TEST(grating_test_471){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 2, 7);
    }

    TEST(grating_test_472){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 3, 1);
    }

    TEST(grating_test_473){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 3, 1);
    }

    TEST(grating_test_474){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 3, 2);
    }

    TEST(grating_test_475){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 3, 2);
    }

    TEST(grating_test_476){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 3, 3);
    }

    TEST(grating_test_477){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 3, 3);
    }

    TEST(grating_test_478){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 3, 5);
    }

    TEST(grating_test_479){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 3, 5);
    }

    TEST(grating_test_480){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 3, 6);
    }

    TEST(grating_test_481){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 3, 6);
    }

    TEST(grating_test_482){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 3, 7);
    }

    TEST(grating_test_483){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 3, 7);
    }

    TEST(grating_test_484){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 5, 1);
    }

    TEST(grating_test_485){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 5, 1);
    }

    TEST(grating_test_486){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 5, 2);
    }

    TEST(grating_test_487){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 5, 2);
    }

    TEST(grating_test_488){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 5, 3);
    }

    TEST(grating_test_489){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 5, 3);
    }

    TEST(grating_test_490){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 5, 5);
    }

    TEST(grating_test_491){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 5, 5);
    }

    TEST(grating_test_492){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 5, 6);
    }

    TEST(grating_test_493){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 5, 6);
    }

    TEST(grating_test_494){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 5, 7);
    }

    TEST(grating_test_495){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 5, 7);
    }

    TEST(grating_test_496){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 6, 1);
    }

    TEST(grating_test_497){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 6, 1);
    }

    TEST(grating_test_498){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 6, 2);
    }

    TEST(grating_test_499){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 6, 2);
    }

    TEST(grating_test_500){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 6, 3);
    }

    TEST(grating_test_501){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 6, 3);
    }

    TEST(grating_test_502){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 6, 5);
    }

    TEST(grating_test_503){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 6, 5);
    }

    TEST(grating_test_504){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 6, 6);
    }

    TEST(grating_test_505){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 6, 6);
    }

    TEST(grating_test_506){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 6, 7);
    }

    TEST(grating_test_507){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 6, 7);
    }

    TEST(grating_test_508){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 7, 1);
    }

    TEST(grating_test_509){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 7, 1);
    }

    TEST(grating_test_510){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 7, 2);
    }

    TEST(grating_test_511){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 7, 2);
    }

    TEST(grating_test_512){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 7, 3);
    }

    TEST(grating_test_513){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 7, 3);
    }

    TEST(grating_test_514){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 7, 5);
    }

    TEST(grating_test_515){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 7, 5);
    }

    TEST(grating_test_516){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 7, 6);
    }

    TEST(grating_test_517){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 7, 6);
    }

    TEST(grating_test_518){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 7, 7);
    }

    TEST(grating_test_519){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 7, 7);
    }

    TEST(grating_test_520){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 1, 1);
    }

    TEST(grating_test_521){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 1, 1);
    }

    TEST(grating_test_522){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 1, 2);
    }

    TEST(grating_test_523){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 1, 2);
    }

    TEST(grating_test_524){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 1, 3);
    }

    TEST(grating_test_525){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 1, 3);
    }

    TEST(grating_test_526){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 1, 5);
    }

    TEST(grating_test_527){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 1, 5);
    }

    TEST(grating_test_528){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 1, 6);
    }

    TEST(grating_test_529){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 1, 6);
    }

    TEST(grating_test_530){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 1, 7);
    }

    TEST(grating_test_531){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 1, 7);
    }

    TEST(grating_test_532){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 2, 1);
    }

    TEST(grating_test_533){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 2, 1);
    }

    TEST(grating_test_534){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 2, 2);
    }

    TEST(grating_test_535){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 2, 2);
    }

    TEST(grating_test_536){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 2, 3);
    }

    TEST(grating_test_537){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 2, 3);
    }

    TEST(grating_test_538){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 2, 5);
    }

    TEST(grating_test_539){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 2, 5);
    }

    TEST(grating_test_540){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 2, 6);
    }

    TEST(grating_test_541){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 2, 6);
    }

    TEST(grating_test_542){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 2, 7);
    }

    TEST(grating_test_543){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 2, 7);
    }

    TEST(grating_test_544){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 3, 1);
    }

    TEST(grating_test_545){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 3, 1);
    }

    TEST(grating_test_546){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 3, 2);
    }

    TEST(grating_test_547){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 3, 2);
    }

    TEST(grating_test_548){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 3, 3);
    }

    TEST(grating_test_549){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 3, 3);
    }

    TEST(grating_test_550){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 3, 5);
    }

    TEST(grating_test_551){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 3, 5);
    }

    TEST(grating_test_552){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 3, 6);
    }

    TEST(grating_test_553){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 3, 6);
    }

    TEST(grating_test_554){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 3, 7);
    }

    TEST(grating_test_555){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 3, 7);
    }

    TEST(grating_test_556){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 5, 1);
    }

    TEST(grating_test_557){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 5, 1);
    }

    TEST(grating_test_558){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 5, 2);
    }

    TEST(grating_test_559){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 5, 2);
    }

    TEST(grating_test_560){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 5, 3);
    }

    TEST(grating_test_561){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 5, 3);
    }

    TEST(grating_test_562){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 5, 5);
    }

    TEST(grating_test_563){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 5, 5);
    }

    TEST(grating_test_564){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 5, 6);
    }

    TEST(grating_test_565){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 5, 6);
    }

    TEST(grating_test_566){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 5, 7);
    }

    TEST(grating_test_567){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 5, 7);
    }

    TEST(grating_test_568){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 6, 1);
    }

    TEST(grating_test_569){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 6, 1);
    }

    TEST(grating_test_570){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 6, 2);
    }

    TEST(grating_test_571){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 6, 2);
    }

    TEST(grating_test_572){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 6, 3);
    }

    TEST(grating_test_573){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 6, 3);
    }

    TEST(grating_test_574){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 6, 5);
    }

    TEST(grating_test_575){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 6, 5);
    }

    TEST(grating_test_576){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 6, 6);
    }

    TEST(grating_test_577){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 6, 6);
    }

    TEST(grating_test_578){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 6, 7);
    }

    TEST(grating_test_579){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 6, 7);
    }

    TEST(grating_test_580){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 7, 1);
    }

    TEST(grating_test_581){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 7, 1);
    }

    TEST(grating_test_582){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 7, 2);
    }

    TEST(grating_test_583){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 7, 2);
    }

    TEST(grating_test_584){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 7, 3);
    }

    TEST(grating_test_585){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 7, 3);
    }

    TEST(grating_test_586){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 7, 5);
    }

    TEST(grating_test_587){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 7, 5);
    }

    TEST(grating_test_588){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 7, 6);
    }

    TEST(grating_test_589){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 7, 6);
    }

    TEST(grating_test_590){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 7, 7);
    }

    TEST(grating_test_591){
         runIntegratorGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 7, 7);
    }

    TEST(grating_test_592){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 1, 1);
    }

    TEST(grating_test_593){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 1, 1);
    }

    TEST(grating_test_594){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 1, 2);
    }

    TEST(grating_test_595){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 1, 2);
    }

    TEST(grating_test_596){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 1, 3);
    }

    TEST(grating_test_597){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 1, 3);
    }

    TEST(grating_test_598){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 1, 5);
    }

    TEST(grating_test_599){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 1, 5);
    }

    TEST(grating_test_600){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 1, 6);
    }

    TEST(grating_test_601){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 1, 6);
    }

    TEST(grating_test_602){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 1, 7);
    }

    TEST(grating_test_603){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 1, 7);
    }

    TEST(grating_test_604){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 2, 1);
    }

    TEST(grating_test_605){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 2, 1);
    }

    TEST(grating_test_606){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 2, 2);
    }

    TEST(grating_test_607){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 2, 2);
    }

    TEST(grating_test_608){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 2, 3);
    }

    TEST(grating_test_609){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 2, 3);
    }

    TEST(grating_test_610){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 2, 5);
    }

    TEST(grating_test_611){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 2, 5);
    }

    TEST(grating_test_612){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 2, 6);
    }

    TEST(grating_test_613){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 2, 6);
    }

    TEST(grating_test_614){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 2, 7);
    }

    TEST(grating_test_615){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 2, 7);
    }

    TEST(grating_test_616){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 3, 1);
    }

    TEST(grating_test_617){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 3, 1);
    }

    TEST(grating_test_618){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 3, 2);
    }

    TEST(grating_test_619){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 3, 2);
    }

    TEST(grating_test_620){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 3, 3);
    }

    TEST(grating_test_621){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 3, 3);
    }

    TEST(grating_test_622){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 3, 5);
    }

    TEST(grating_test_623){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 3, 5);
    }

    TEST(grating_test_624){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 3, 6);
    }

    TEST(grating_test_625){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 3, 6);
    }

    TEST(grating_test_626){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 3, 7);
    }

    TEST(grating_test_627){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 3, 7);
    }

    TEST(grating_test_628){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 5, 1);
    }

    TEST(grating_test_629){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 5, 1);
    }

    TEST(grating_test_630){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 5, 2);
    }

    TEST(grating_test_631){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 5, 2);
    }

    TEST(grating_test_632){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 5, 3);
    }

    TEST(grating_test_633){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 5, 3);
    }

    TEST(grating_test_634){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 5, 5);
    }

    TEST(grating_test_635){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 5, 5);
    }

    TEST(grating_test_636){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 5, 6);
    }

    TEST(grating_test_637){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 5, 6);
    }

    TEST(grating_test_638){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 5, 7);
    }

    TEST(grating_test_639){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 5, 7);
    }

    TEST(grating_test_640){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 6, 1);
    }

    TEST(grating_test_641){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 6, 1);
    }

    TEST(grating_test_642){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 6, 2);
    }

    TEST(grating_test_643){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 6, 2);
    }

    TEST(grating_test_644){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 6, 3);
    }

    TEST(grating_test_645){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 6, 3);
    }

    TEST(grating_test_646){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 6, 5);
    }

    TEST(grating_test_647){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 6, 5);
    }

    TEST(grating_test_648){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 6, 6);
    }

    TEST(grating_test_649){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 6, 6);
    }

    TEST(grating_test_650){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 6, 7);
    }

    TEST(grating_test_651){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 6, 7);
    }

    TEST(grating_test_652){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 7, 1);
    }

    TEST(grating_test_653){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 7, 1);
    }

    TEST(grating_test_654){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 7, 2);
    }

    TEST(grating_test_655){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 7, 2);
    }

    TEST(grating_test_656){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 7, 3);
    }

    TEST(grating_test_657){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 7, 3);
    }

    TEST(grating_test_658){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 7, 5);
    }

    TEST(grating_test_659){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 7, 5);
    }

    TEST(grating_test_660){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 7, 6);
    }

    TEST(grating_test_661){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 7, 6);
    }

    TEST(grating_test_662){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 7, 7);
    }

    TEST(grating_test_663){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 7, 7);
    }

    TEST(grating_test_664){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 1, 1);
    }

    TEST(grating_test_665){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 1, 1);
    }

    TEST(grating_test_666){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 1, 2);
    }

    TEST(grating_test_667){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 1, 2);
    }

    TEST(grating_test_668){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 1, 3);
    }

    TEST(grating_test_669){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 1, 3);
    }

    TEST(grating_test_670){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 1, 5);
    }

    TEST(grating_test_671){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 1, 5);
    }

    TEST(grating_test_672){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 1, 6);
    }

    TEST(grating_test_673){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 1, 6);
    }

    TEST(grating_test_674){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 1, 7);
    }

    TEST(grating_test_675){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 1, 7);
    }

    TEST(grating_test_676){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 2, 1);
    }

    TEST(grating_test_677){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 2, 1);
    }

    TEST(grating_test_678){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 2, 2);
    }

    TEST(grating_test_679){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 2, 2);
    }

    TEST(grating_test_680){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 2, 3);
    }

    TEST(grating_test_681){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 2, 3);
    }

    TEST(grating_test_682){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 2, 5);
    }

    TEST(grating_test_683){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 2, 5);
    }

    TEST(grating_test_684){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 2, 6);
    }

    TEST(grating_test_685){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 2, 6);
    }

    TEST(grating_test_686){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 2, 7);
    }

    TEST(grating_test_687){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 2, 7);
    }

    TEST(grating_test_688){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 3, 1);
    }

    TEST(grating_test_689){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 3, 1);
    }

    TEST(grating_test_690){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 3, 2);
    }

    TEST(grating_test_691){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 3, 2);
    }

    TEST(grating_test_692){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 3, 3);
    }

    TEST(grating_test_693){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 3, 3);
    }

    TEST(grating_test_694){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 3, 5);
    }

    TEST(grating_test_695){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 3, 5);
    }

    TEST(grating_test_696){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 3, 6);
    }

    TEST(grating_test_697){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 3, 6);
    }

    TEST(grating_test_698){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 3, 7);
    }

    TEST(grating_test_699){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 3, 7);
    }

    TEST(grating_test_700){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 5, 1);
    }

    TEST(grating_test_701){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 5, 1);
    }

    TEST(grating_test_702){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 5, 2);
    }

    TEST(grating_test_703){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 5, 2);
    }

    TEST(grating_test_704){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 5, 3);
    }

    TEST(grating_test_705){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 5, 3);
    }

    TEST(grating_test_706){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 5, 5);
    }

    TEST(grating_test_707){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 5, 5);
    }

    TEST(grating_test_708){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 5, 6);
    }

    TEST(grating_test_709){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 5, 6);
    }

    TEST(grating_test_710){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 5, 7);
    }

    TEST(grating_test_711){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 5, 7);
    }

    TEST(grating_test_712){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 6, 1);
    }

    TEST(grating_test_713){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 6, 1);
    }

    TEST(grating_test_714){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 6, 2);
    }

    TEST(grating_test_715){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 6, 2);
    }

    TEST(grating_test_716){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 6, 3);
    }

    TEST(grating_test_717){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 6, 3);
    }

    TEST(grating_test_718){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 6, 5);
    }

    TEST(grating_test_719){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 6, 5);
    }

    TEST(grating_test_720){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 6, 6);
    }

    TEST(grating_test_721){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 6, 6);
    }

    TEST(grating_test_722){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 6, 7);
    }

    TEST(grating_test_723){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 6, 7);
    }

    TEST(grating_test_724){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 7, 1);
    }

    TEST(grating_test_725){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 7, 1);
    }

    TEST(grating_test_726){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 7, 2);
    }

    TEST(grating_test_727){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 7, 2);
    }

    TEST(grating_test_728){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 7, 3);
    }

    TEST(grating_test_729){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 7, 3);
    }

    TEST(grating_test_730){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 7, 5);
    }

    TEST(grating_test_731){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 7, 5);
    }

    TEST(grating_test_732){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 7, 6);
    }

    TEST(grating_test_733){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 7, 6);
    }

    TEST(grating_test_734){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 7, 7);
    }

    TEST(grating_test_735){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 7, 7);
    }

    TEST(grating_test_736){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 1, 1);
    }

    TEST(grating_test_737){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 1, 1);
    }

    TEST(grating_test_738){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 1, 2);
    }

    TEST(grating_test_739){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 1, 2);
    }

    TEST(grating_test_740){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 1, 3);
    }

    TEST(grating_test_741){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 1, 3);
    }

    TEST(grating_test_742){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 1, 5);
    }

    TEST(grating_test_743){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 1, 5);
    }

    TEST(grating_test_744){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 1, 6);
    }

    TEST(grating_test_745){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 1, 6);
    }

    TEST(grating_test_746){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 1, 7);
    }

    TEST(grating_test_747){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 1, 7);
    }

    TEST(grating_test_748){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 2, 1);
    }

    TEST(grating_test_749){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 2, 1);
    }

    TEST(grating_test_750){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 2, 2);
    }

    TEST(grating_test_751){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 2, 2);
    }

    TEST(grating_test_752){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 2, 3);
    }

    TEST(grating_test_753){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 2, 3);
    }

    TEST(grating_test_754){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 2, 5);
    }

    TEST(grating_test_755){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 2, 5);
    }

    TEST(grating_test_756){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 2, 6);
    }

    TEST(grating_test_757){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 2, 6);
    }

    TEST(grating_test_758){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 2, 7);
    }

    TEST(grating_test_759){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 2, 7);
    }

    TEST(grating_test_760){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 3, 1);
    }

    TEST(grating_test_761){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 3, 1);
    }

    TEST(grating_test_762){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 3, 2);
    }

    TEST(grating_test_763){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 3, 2);
    }

    TEST(grating_test_764){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 3, 3);
    }

    TEST(grating_test_765){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 3, 3);
    }

    TEST(grating_test_766){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 3, 5);
    }

    TEST(grating_test_767){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 3, 5);
    }

    TEST(grating_test_768){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 3, 6);
    }

    TEST(grating_test_769){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 3, 6);
    }

    TEST(grating_test_770){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 3, 7);
    }

    TEST(grating_test_771){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 3, 7);
    }

    TEST(grating_test_772){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 5, 1);
    }

    TEST(grating_test_773){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 5, 1);
    }

    TEST(grating_test_774){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 5, 2);
    }

    TEST(grating_test_775){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 5, 2);
    }

    TEST(grating_test_776){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 5, 3);
    }

    TEST(grating_test_777){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 5, 3);
    }

    TEST(grating_test_778){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 5, 5);
    }

    TEST(grating_test_779){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 5, 5);
    }

    TEST(grating_test_780){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 5, 6);
    }

    TEST(grating_test_781){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 5, 6);
    }

    TEST(grating_test_782){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 5, 7);
    }

    TEST(grating_test_783){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 5, 7);
    }

    TEST(grating_test_784){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 6, 1);
    }

    TEST(grating_test_785){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 6, 1);
    }

    TEST(grating_test_786){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 6, 2);
    }

    TEST(grating_test_787){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 6, 2);
    }

    TEST(grating_test_788){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 6, 3);
    }

    TEST(grating_test_789){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 6, 3);
    }

    TEST(grating_test_790){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 6, 5);
    }

    TEST(grating_test_791){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 6, 5);
    }

    TEST(grating_test_792){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 6, 6);
    }

    TEST(grating_test_793){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 6, 6);
    }

    TEST(grating_test_794){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 6, 7);
    }

    TEST(grating_test_795){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 6, 7);
    }

    TEST(grating_test_796){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 7, 1);
    }

    TEST(grating_test_797){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 7, 1);
    }

    TEST(grating_test_798){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 7, 2);
    }

    TEST(grating_test_799){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 7, 2);
    }

    TEST(grating_test_800){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 7, 3);
    }

    TEST(grating_test_801){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 7, 3);
    }

    TEST(grating_test_802){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 7, 5);
    }

    TEST(grating_test_803){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 7, 5);
    }

    TEST(grating_test_804){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 7, 6);
    }

    TEST(grating_test_805){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 7, 6);
    }

    TEST(grating_test_806){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 7, 7);
    }

    TEST(grating_test_807){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 7, 7);
    }

    TEST(grating_test_808){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 1, 1);
    }

    TEST(grating_test_809){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 1, 1);
    }

    TEST(grating_test_810){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 1, 2);
    }

    TEST(grating_test_811){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 1, 2);
    }

    TEST(grating_test_812){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 1, 3);
    }

    TEST(grating_test_813){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 1, 3);
    }

    TEST(grating_test_814){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 1, 5);
    }

    TEST(grating_test_815){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 1, 5);
    }

    TEST(grating_test_816){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 1, 6);
    }

    TEST(grating_test_817){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 1, 6);
    }

    TEST(grating_test_818){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 1, 7);
    }

    TEST(grating_test_819){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 1, 7);
    }

    TEST(grating_test_820){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 2, 1);
    }

    TEST(grating_test_821){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 2, 1);
    }

    TEST(grating_test_822){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 2, 2);
    }

    TEST(grating_test_823){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 2, 2);
    }

    TEST(grating_test_824){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 2, 3);
    }

    TEST(grating_test_825){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 2, 3);
    }

    TEST(grating_test_826){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 2, 5);
    }

    TEST(grating_test_827){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 2, 5);
    }

    TEST(grating_test_828){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 2, 6);
    }

    TEST(grating_test_829){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 2, 6);
    }

    TEST(grating_test_830){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 2, 7);
    }

    TEST(grating_test_831){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 2, 7);
    }

    TEST(grating_test_832){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 3, 1);
    }

    TEST(grating_test_833){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 3, 1);
    }

    TEST(grating_test_834){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 3, 2);
    }

    TEST(grating_test_835){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 3, 2);
    }

    TEST(grating_test_836){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 3, 3);
    }

    TEST(grating_test_837){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 3, 3);
    }

    TEST(grating_test_838){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 3, 5);
    }

    TEST(grating_test_839){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 3, 5);
    }

    TEST(grating_test_840){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 3, 6);
    }

    TEST(grating_test_841){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 3, 6);
    }

    TEST(grating_test_842){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 3, 7);
    }

    TEST(grating_test_843){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 3, 7);
    }

    TEST(grating_test_844){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 5, 1);
    }

    TEST(grating_test_845){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 5, 1);
    }

    TEST(grating_test_846){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 5, 2);
    }

    TEST(grating_test_847){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 5, 2);
    }

    TEST(grating_test_848){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 5, 3);
    }

    TEST(grating_test_849){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 5, 3);
    }

    TEST(grating_test_850){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 5, 5);
    }

    TEST(grating_test_851){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 5, 5);
    }

    TEST(grating_test_852){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 5, 6);
    }

    TEST(grating_test_853){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 5, 6);
    }

    TEST(grating_test_854){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 5, 7);
    }

    TEST(grating_test_855){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 5, 7);
    }

    TEST(grating_test_856){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 6, 1);
    }

    TEST(grating_test_857){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 6, 1);
    }

    TEST(grating_test_858){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 6, 2);
    }

    TEST(grating_test_859){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 6, 2);
    }

    TEST(grating_test_860){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 6, 3);
    }

    TEST(grating_test_861){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 6, 3);
    }

    TEST(grating_test_862){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 6, 5);
    }

    TEST(grating_test_863){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 6, 5);
    }

    TEST(grating_test_864){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 6, 6);
    }

    TEST(grating_test_865){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 6, 6);
    }

    TEST(grating_test_866){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 6, 7);
    }

    TEST(grating_test_867){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 6, 7);
    }

    TEST(grating_test_868){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 7, 1);
    }

    TEST(grating_test_869){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 7, 1);
    }

    TEST(grating_test_870){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 7, 2);
    }

    TEST(grating_test_871){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 7, 2);
    }

    TEST(grating_test_872){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 7, 3);
    }

    TEST(grating_test_873){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 7, 3);
    }

    TEST(grating_test_874){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 7, 5);
    }

    TEST(grating_test_875){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 7, 5);
    }

    TEST(grating_test_876){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 7, 6);
    }

    TEST(grating_test_877){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 7, 6);
    }

    TEST(grating_test_878){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 7, 7);
    }

    TEST(grating_test_879){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 7, 7);
    }

    TEST(grating_test_880){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 1, 1);
    }

    TEST(grating_test_881){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 1, 1);
    }

    TEST(grating_test_882){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 1, 2);
    }

    TEST(grating_test_883){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 1, 2);
    }

    TEST(grating_test_884){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 1, 3);
    }

    TEST(grating_test_885){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 1, 3);
    }

    TEST(grating_test_886){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 1, 5);
    }

    TEST(grating_test_887){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 1, 5);
    }

    TEST(grating_test_888){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 1, 6);
    }

    TEST(grating_test_889){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 1, 6);
    }

    TEST(grating_test_890){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 1, 7);
    }

    TEST(grating_test_891){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 1, 7);
    }

    TEST(grating_test_892){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 2, 1);
    }

    TEST(grating_test_893){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 2, 1);
    }

    TEST(grating_test_894){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 2, 2);
    }

    TEST(grating_test_895){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 2, 2);
    }

    TEST(grating_test_896){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 2, 3);
    }

    TEST(grating_test_897){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 2, 3);
    }

    TEST(grating_test_898){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 2, 5);
    }

    TEST(grating_test_899){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 2, 5);
    }

    TEST(grating_test_900){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 2, 6);
    }

    TEST(grating_test_901){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 2, 6);
    }

    TEST(grating_test_902){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 2, 7);
    }

    TEST(grating_test_903){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 2, 7);
    }

    TEST(grating_test_904){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 3, 1);
    }

    TEST(grating_test_905){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 3, 1);
    }

    TEST(grating_test_906){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 3, 2);
    }

    TEST(grating_test_907){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 3, 2);
    }

    TEST(grating_test_908){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 3, 3);
    }

    TEST(grating_test_909){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 3, 3);
    }

    TEST(grating_test_910){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 3, 5);
    }

    TEST(grating_test_911){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 3, 5);
    }

    TEST(grating_test_912){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 3, 6);
    }

    TEST(grating_test_913){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 3, 6);
    }

    TEST(grating_test_914){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 3, 7);
    }

    TEST(grating_test_915){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 3, 7);
    }

    TEST(grating_test_916){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 5, 1);
    }

    TEST(grating_test_917){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 5, 1);
    }

    TEST(grating_test_918){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 5, 2);
    }

    TEST(grating_test_919){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 5, 2);
    }

    TEST(grating_test_920){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 5, 3);
    }

    TEST(grating_test_921){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 5, 3);
    }

    TEST(grating_test_922){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 5, 5);
    }

    TEST(grating_test_923){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 5, 5);
    }

    TEST(grating_test_924){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 5, 6);
    }

    TEST(grating_test_925){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 5, 6);
    }

    TEST(grating_test_926){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 5, 7);
    }

    TEST(grating_test_927){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 5, 7);
    }

    TEST(grating_test_928){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 6, 1);
    }

    TEST(grating_test_929){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 6, 1);
    }

    TEST(grating_test_930){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 6, 2);
    }

    TEST(grating_test_931){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 6, 2);
    }

    TEST(grating_test_932){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 6, 3);
    }

    TEST(grating_test_933){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 6, 3);
    }

    TEST(grating_test_934){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 6, 5);
    }

    TEST(grating_test_935){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 6, 5);
    }

    TEST(grating_test_936){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 6, 6);
    }

    TEST(grating_test_937){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 6, 6);
    }

    TEST(grating_test_938){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 6, 7);
    }

    TEST(grating_test_939){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 6, 7);
    }

    TEST(grating_test_940){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 7, 1);
    }

    TEST(grating_test_941){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 7, 1);
    }

    TEST(grating_test_942){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 7, 2);
    }

    TEST(grating_test_943){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 7, 2);
    }

    TEST(grating_test_944){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 7, 3);
    }

    TEST(grating_test_945){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 7, 3);
    }

    TEST(grating_test_946){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 7, 5);
    }

    TEST(grating_test_947){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 7, 5);
    }

    TEST(grating_test_948){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 7, 6);
    }

    TEST(grating_test_949){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 7, 6);
    }

    TEST(grating_test_950){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 7, 7);
    }

    TEST(grating_test_951){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 7, 7);
    }

    TEST(grating_test_952){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 1, 1);
    }

    TEST(grating_test_953){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 1, 1);
    }

    TEST(grating_test_954){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 1, 2);
    }

    TEST(grating_test_955){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 1, 2);
    }

    TEST(grating_test_956){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 1, 3);
    }

    TEST(grating_test_957){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 1, 3);
    }

    TEST(grating_test_958){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 1, 5);
    }

    TEST(grating_test_959){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 1, 5);
    }

    TEST(grating_test_960){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 1, 6);
    }

    TEST(grating_test_961){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 1, 6);
    }

    TEST(grating_test_962){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 1, 7);
    }

    TEST(grating_test_963){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 1, 7);
    }

    TEST(grating_test_964){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 2, 1);
    }

    TEST(grating_test_965){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 2, 1);
    }

    TEST(grating_test_966){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 2, 2);
    }

    TEST(grating_test_967){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 2, 2);
    }

    TEST(grating_test_968){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 2, 3);
    }

    TEST(grating_test_969){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 2, 3);
    }

    TEST(grating_test_970){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 2, 5);
    }

    TEST(grating_test_971){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 2, 5);
    }

    TEST(grating_test_972){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 2, 6);
    }

    TEST(grating_test_973){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 2, 6);
    }

    TEST(grating_test_974){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 2, 7);
    }

    TEST(grating_test_975){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 2, 7);
    }

    TEST(grating_test_976){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 3, 1);
    }

    TEST(grating_test_977){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 3, 1);
    }

    TEST(grating_test_978){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 3, 2);
    }

    TEST(grating_test_979){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 3, 2);
    }

    TEST(grating_test_980){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 3, 3);
    }

    TEST(grating_test_981){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 3, 3);
    }

    TEST(grating_test_982){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 3, 5);
    }

    TEST(grating_test_983){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 3, 5);
    }

    TEST(grating_test_984){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 3, 6);
    }

    TEST(grating_test_985){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 3, 6);
    }

    TEST(grating_test_986){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 3, 7);
    }

    TEST(grating_test_987){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 3, 7);
    }

    TEST(grating_test_988){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 5, 1);
    }

    TEST(grating_test_989){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 5, 1);
    }

    TEST(grating_test_990){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 5, 2);
    }

    TEST(grating_test_991){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 5, 2);
    }

    TEST(grating_test_992){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 5, 3);
    }

    TEST(grating_test_993){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 5, 3);
    }

    TEST(grating_test_994){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 5, 5);
    }

    TEST(grating_test_995){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 5, 5);
    }

    TEST(grating_test_996){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 5, 6);
    }

    TEST(grating_test_997){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 5, 6);
    }

    TEST(grating_test_998){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 5, 7);
    }

    TEST(grating_test_999){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 5, 7);
    }

    TEST(grating_test_1000){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 6, 1);
    }

    TEST(grating_test_1001){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 6, 1);
    }

    TEST(grating_test_1002){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 6, 2);
    }

    TEST(grating_test_1003){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 6, 2);
    }

    TEST(grating_test_1004){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 6, 3);
    }

    TEST(grating_test_1005){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 6, 3);
    }

    TEST(grating_test_1006){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 6, 5);
    }

    TEST(grating_test_1007){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 6, 5);
    }

    TEST(grating_test_1008){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 6, 6);
    }

    TEST(grating_test_1009){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 6, 6);
    }

    TEST(grating_test_1010){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 6, 7);
    }

    TEST(grating_test_1011){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 6, 7);
    }

    TEST(grating_test_1012){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 7, 1);
    }

    TEST(grating_test_1013){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 7, 1);
    }

    TEST(grating_test_1014){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 7, 2);
    }

    TEST(grating_test_1015){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 7, 2);
    }

    TEST(grating_test_1016){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 7, 3);
    }

    TEST(grating_test_1017){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 7, 3);
    }

    TEST(grating_test_1018){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 7, 5);
    }

    TEST(grating_test_1019){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 7, 5);
    }

    TEST(grating_test_1020){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 7, 6);
    }

    TEST(grating_test_1021){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 7, 6);
    }

    TEST(grating_test_1022){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 7, 7);
    }

    TEST(grating_test_1023){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 7, 7);
    }

    TEST(grating_test_1024){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 1, 1);
    }

    TEST(grating_test_1025){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 1, 1);
    }

    TEST(grating_test_1026){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 1, 2);
    }

    TEST(grating_test_1027){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 1, 2);
    }

    TEST(grating_test_1028){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 1, 3);
    }

    TEST(grating_test_1029){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 1, 3);
    }

    TEST(grating_test_1030){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 1, 5);
    }

    TEST(grating_test_1031){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 1, 5);
    }

    TEST(grating_test_1032){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 1, 6);
    }

    TEST(grating_test_1033){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 1, 6);
    }

    TEST(grating_test_1034){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 1, 7);
    }

    TEST(grating_test_1035){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 1, 7);
    }

    TEST(grating_test_1036){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 2, 1);
    }

    TEST(grating_test_1037){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 2, 1);
    }

    TEST(grating_test_1038){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 2, 2);
    }

    TEST(grating_test_1039){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 2, 2);
    }

    TEST(grating_test_1040){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 2, 3);
    }

    TEST(grating_test_1041){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 2, 3);
    }

    TEST(grating_test_1042){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 2, 5);
    }

    TEST(grating_test_1043){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 2, 5);
    }

    TEST(grating_test_1044){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 2, 6);
    }

    TEST(grating_test_1045){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 2, 6);
    }

    TEST(grating_test_1046){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 2, 7);
    }

    TEST(grating_test_1047){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 2, 7);
    }

    TEST(grating_test_1048){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 3, 1);
    }

    TEST(grating_test_1049){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 3, 1);
    }

    TEST(grating_test_1050){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 3, 2);
    }

    TEST(grating_test_1051){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 3, 2);
    }

    TEST(grating_test_1052){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 3, 3);
    }

    TEST(grating_test_1053){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 3, 3);
    }

    TEST(grating_test_1054){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 3, 5);
    }

    TEST(grating_test_1055){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 3, 5);
    }

    TEST(grating_test_1056){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 3, 6);
    }

    TEST(grating_test_1057){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 3, 6);
    }

    TEST(grating_test_1058){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 3, 7);
    }

    TEST(grating_test_1059){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 3, 7);
    }

    TEST(grating_test_1060){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 5, 1);
    }

    TEST(grating_test_1061){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 5, 1);
    }

    TEST(grating_test_1062){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 5, 2);
    }

    TEST(grating_test_1063){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 5, 2);
    }

    TEST(grating_test_1064){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 5, 3);
    }

    TEST(grating_test_1065){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 5, 3);
    }

    TEST(grating_test_1066){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 5, 5);
    }

    TEST(grating_test_1067){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 5, 5);
    }

    TEST(grating_test_1068){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 5, 6);
    }

    TEST(grating_test_1069){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 5, 6);
    }

    TEST(grating_test_1070){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 5, 7);
    }

    TEST(grating_test_1071){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 5, 7);
    }

    TEST(grating_test_1072){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 6, 1);
    }

    TEST(grating_test_1073){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 6, 1);
    }

    TEST(grating_test_1074){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 6, 2);
    }

    TEST(grating_test_1075){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 6, 2);
    }

    TEST(grating_test_1076){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 6, 3);
    }

    TEST(grating_test_1077){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 6, 3);
    }

    TEST(grating_test_1078){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 6, 5);
    }

    TEST(grating_test_1079){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 6, 5);
    }

    TEST(grating_test_1080){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 6, 6);
    }

    TEST(grating_test_1081){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 6, 6);
    }

    TEST(grating_test_1082){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 6, 7);
    }

    TEST(grating_test_1083){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 6, 7);
    }

    TEST(grating_test_1084){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 7, 1);
    }

    TEST(grating_test_1085){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 7, 1);
    }

    TEST(grating_test_1086){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 7, 2);
    }

    TEST(grating_test_1087){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 7, 2);
    }

    TEST(grating_test_1088){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 7, 3);
    }

    TEST(grating_test_1089){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 7, 3);
    }

    TEST(grating_test_1090){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 7, 5);
    }

    TEST(grating_test_1091){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 7, 5);
    }

    TEST(grating_test_1092){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 7, 6);
    }

    TEST(grating_test_1093){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 7, 6);
    }

    TEST(grating_test_1094){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 7, 7);
    }

    TEST(grating_test_1095){
         runIntegratorGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 7, 7);
    }

    TEST(grating_test_1096){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 1, 1);
    }

    TEST(grating_test_1097){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 1, 1);
    }

    TEST(grating_test_1098){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 1, 2);
    }

    TEST(grating_test_1099){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 1, 2);
    }

    TEST(grating_test_1100){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 1, 3);
    }

    TEST(grating_test_1101){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 1, 3);
    }

    TEST(grating_test_1102){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 1, 5);
    }

    TEST(grating_test_1103){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 1, 5);
    }

    TEST(grating_test_1104){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 1, 6);
    }

    TEST(grating_test_1105){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 1, 6);
    }

    TEST(grating_test_1106){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 1, 7);
    }

    TEST(grating_test_1107){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 1, 7);
    }

    TEST(grating_test_1108){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 2, 1);
    }

    TEST(grating_test_1109){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 2, 1);
    }

    TEST(grating_test_1110){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 2, 2);
    }

    TEST(grating_test_1111){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 2, 2);
    }

    TEST(grating_test_1112){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 2, 3);
    }

    TEST(grating_test_1113){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 2, 3);
    }

    TEST(grating_test_1114){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 2, 5);
    }

    TEST(grating_test_1115){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 2, 5);
    }

    TEST(grating_test_1116){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 2, 6);
    }

    TEST(grating_test_1117){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 2, 6);
    }

    TEST(grating_test_1118){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 2, 7);
    }

    TEST(grating_test_1119){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 2, 7);
    }

    TEST(grating_test_1120){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 3, 1);
    }

    TEST(grating_test_1121){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 3, 1);
    }

    TEST(grating_test_1122){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 3, 2);
    }

    TEST(grating_test_1123){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 3, 2);
    }

    TEST(grating_test_1124){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 3, 3);
    }

    TEST(grating_test_1125){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 3, 3);
    }

    TEST(grating_test_1126){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 3, 5);
    }

    TEST(grating_test_1127){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 3, 5);
    }

    TEST(grating_test_1128){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 3, 6);
    }

    TEST(grating_test_1129){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 3, 6);
    }

    TEST(grating_test_1130){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 3, 7);
    }

    TEST(grating_test_1131){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 3, 7);
    }

    TEST(grating_test_1132){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 5, 1);
    }

    TEST(grating_test_1133){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 5, 1);
    }

    TEST(grating_test_1134){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 5, 2);
    }

    TEST(grating_test_1135){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 5, 2);
    }

    TEST(grating_test_1136){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 5, 3);
    }

    TEST(grating_test_1137){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 5, 3);
    }

    TEST(grating_test_1138){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 5, 5);
    }

    TEST(grating_test_1139){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 5, 5);
    }

    TEST(grating_test_1140){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 5, 6);
    }

    TEST(grating_test_1141){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 5, 6);
    }

    TEST(grating_test_1142){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 5, 7);
    }

    TEST(grating_test_1143){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 5, 7);
    }

    TEST(grating_test_1144){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 6, 1);
    }

    TEST(grating_test_1145){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 6, 1);
    }

    TEST(grating_test_1146){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 6, 2);
    }

    TEST(grating_test_1147){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 6, 2);
    }

    TEST(grating_test_1148){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 6, 3);
    }

    TEST(grating_test_1149){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 6, 3);
    }

    TEST(grating_test_1150){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 6, 5);
    }

    TEST(grating_test_1151){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 6, 5);
    }

    TEST(grating_test_1152){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 6, 6);
    }

    TEST(grating_test_1153){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 6, 6);
    }

    TEST(grating_test_1154){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 6, 7);
    }

    TEST(grating_test_1155){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 6, 7);
    }

    TEST(grating_test_1156){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 7, 1);
    }

    TEST(grating_test_1157){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 7, 1);
    }

    TEST(grating_test_1158){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 7, 2);
    }

    TEST(grating_test_1159){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 7, 2);
    }

    TEST(grating_test_1160){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 7, 3);
    }

    TEST(grating_test_1161){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 7, 3);
    }

    TEST(grating_test_1162){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 7, 5);
    }

    TEST(grating_test_1163){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 7, 5);
    }

    TEST(grating_test_1164){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 7, 6);
    }

    TEST(grating_test_1165){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 7, 6);
    }

    TEST(grating_test_1166){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 7, 7);
    }

    TEST(grating_test_1167){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 7, 7);
    }

    TEST(grating_test_1168){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 1, 1);
    }

    TEST(grating_test_1169){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 1, 1);
    }

    TEST(grating_test_1170){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 1, 2);
    }

    TEST(grating_test_1171){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 1, 2);
    }

    TEST(grating_test_1172){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 1, 3);
    }

    TEST(grating_test_1173){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 1, 3);
    }

    TEST(grating_test_1174){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 1, 5);
    }

    TEST(grating_test_1175){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 1, 5);
    }

    TEST(grating_test_1176){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 1, 6);
    }

    TEST(grating_test_1177){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 1, 6);
    }

    TEST(grating_test_1178){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 1, 7);
    }

    TEST(grating_test_1179){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 1, 7);
    }

    TEST(grating_test_1180){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 2, 1);
    }

    TEST(grating_test_1181){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 2, 1);
    }

    TEST(grating_test_1182){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 2, 2);
    }

    TEST(grating_test_1183){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 2, 2);
    }

    TEST(grating_test_1184){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 2, 3);
    }

    TEST(grating_test_1185){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 2, 3);
    }

    TEST(grating_test_1186){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 2, 5);
    }

    TEST(grating_test_1187){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 2, 5);
    }

    TEST(grating_test_1188){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 2, 6);
    }

    TEST(grating_test_1189){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 2, 6);
    }

    TEST(grating_test_1190){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 2, 7);
    }

    TEST(grating_test_1191){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 2, 7);
    }

    TEST(grating_test_1192){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 3, 1);
    }

    TEST(grating_test_1193){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 3, 1);
    }

    TEST(grating_test_1194){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 3, 2);
    }

    TEST(grating_test_1195){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 3, 2);
    }

    TEST(grating_test_1196){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 3, 3);
    }

    TEST(grating_test_1197){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 3, 3);
    }

    TEST(grating_test_1198){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 3, 5);
    }

    TEST(grating_test_1199){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 3, 5);
    }

    TEST(grating_test_1200){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 3, 6);
    }

    TEST(grating_test_1201){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 3, 6);
    }

    TEST(grating_test_1202){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 3, 7);
    }

    TEST(grating_test_1203){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 3, 7);
    }

    TEST(grating_test_1204){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 5, 1);
    }

    TEST(grating_test_1205){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 5, 1);
    }

    TEST(grating_test_1206){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 5, 2);
    }

    TEST(grating_test_1207){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 5, 2);
    }

    TEST(grating_test_1208){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 5, 3);
    }

    TEST(grating_test_1209){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 5, 3);
    }

    TEST(grating_test_1210){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 5, 5);
    }

    TEST(grating_test_1211){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 5, 5);
    }

    TEST(grating_test_1212){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 5, 6);
    }

    TEST(grating_test_1213){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 5, 6);
    }

    TEST(grating_test_1214){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 5, 7);
    }

    TEST(grating_test_1215){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 5, 7);
    }

    TEST(grating_test_1216){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 6, 1);
    }

    TEST(grating_test_1217){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 6, 1);
    }

    TEST(grating_test_1218){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 6, 2);
    }

    TEST(grating_test_1219){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 6, 2);
    }

    TEST(grating_test_1220){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 6, 3);
    }

    TEST(grating_test_1221){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 6, 3);
    }

    TEST(grating_test_1222){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 6, 5);
    }

    TEST(grating_test_1223){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 6, 5);
    }

    TEST(grating_test_1224){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 6, 6);
    }

    TEST(grating_test_1225){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 6, 6);
    }

    TEST(grating_test_1226){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 6, 7);
    }

    TEST(grating_test_1227){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 6, 7);
    }

    TEST(grating_test_1228){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 7, 1);
    }

    TEST(grating_test_1229){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 7, 1);
    }

    TEST(grating_test_1230){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 7, 2);
    }

    TEST(grating_test_1231){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 7, 2);
    }

    TEST(grating_test_1232){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 7, 3);
    }

    TEST(grating_test_1233){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 7, 3);
    }

    TEST(grating_test_1234){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 7, 5);
    }

    TEST(grating_test_1235){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 7, 5);
    }

    TEST(grating_test_1236){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 7, 6);
    }

    TEST(grating_test_1237){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 7, 6);
    }

    TEST(grating_test_1238){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 7, 7);
    }

    TEST(grating_test_1239){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 7, 7);
    }

    TEST(grating_test_1240){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 1, 1);
    }

    TEST(grating_test_1241){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 1, 1);
    }

    TEST(grating_test_1242){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 1, 2);
    }

    TEST(grating_test_1243){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 1, 2);
    }

    TEST(grating_test_1244){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 1, 3);
    }

    TEST(grating_test_1245){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 1, 3);
    }

    TEST(grating_test_1246){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 1, 5);
    }

    TEST(grating_test_1247){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 1, 5);
    }

    TEST(grating_test_1248){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 1, 6);
    }

    TEST(grating_test_1249){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 1, 6);
    }

    TEST(grating_test_1250){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 1, 7);
    }

    TEST(grating_test_1251){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 1, 7);
    }

    TEST(grating_test_1252){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 2, 1);
    }

    TEST(grating_test_1253){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 2, 1);
    }

    TEST(grating_test_1254){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 2, 2);
    }

    TEST(grating_test_1255){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 2, 2);
    }

    TEST(grating_test_1256){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 2, 3);
    }

    TEST(grating_test_1257){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 2, 3);
    }

    TEST(grating_test_1258){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 2, 5);
    }

    TEST(grating_test_1259){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 2, 5);
    }

    TEST(grating_test_1260){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 2, 6);
    }

    TEST(grating_test_1261){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 2, 6);
    }

    TEST(grating_test_1262){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 2, 7);
    }

    TEST(grating_test_1263){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 2, 7);
    }

    TEST(grating_test_1264){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 3, 1);
    }

    TEST(grating_test_1265){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 3, 1);
    }

    TEST(grating_test_1266){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 3, 2);
    }

    TEST(grating_test_1267){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 3, 2);
    }

    TEST(grating_test_1268){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 3, 3);
    }

    TEST(grating_test_1269){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 3, 3);
    }

    TEST(grating_test_1270){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 3, 5);
    }

    TEST(grating_test_1271){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 3, 5);
    }

    TEST(grating_test_1272){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 3, 6);
    }

    TEST(grating_test_1273){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 3, 6);
    }

    TEST(grating_test_1274){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 3, 7);
    }

    TEST(grating_test_1275){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 3, 7);
    }

    TEST(grating_test_1276){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 5, 1);
    }

    TEST(grating_test_1277){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 5, 1);
    }

    TEST(grating_test_1278){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 5, 2);
    }

    TEST(grating_test_1279){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 5, 2);
    }

    TEST(grating_test_1280){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 5, 3);
    }

    TEST(grating_test_1281){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 5, 3);
    }

    TEST(grating_test_1282){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 5, 5);
    }

    TEST(grating_test_1283){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 5, 5);
    }

    TEST(grating_test_1284){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 5, 6);
    }

    TEST(grating_test_1285){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 5, 6);
    }

    TEST(grating_test_1286){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 5, 7);
    }

    TEST(grating_test_1287){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 5, 7);
    }

    TEST(grating_test_1288){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 6, 1);
    }

    TEST(grating_test_1289){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 6, 1);
    }

    TEST(grating_test_1290){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 6, 2);
    }

    TEST(grating_test_1291){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 6, 2);
    }

    TEST(grating_test_1292){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 6, 3);
    }

    TEST(grating_test_1293){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 6, 3);
    }

    TEST(grating_test_1294){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 6, 5);
    }

    TEST(grating_test_1295){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 6, 5);
    }

    TEST(grating_test_1296){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 6, 6);
    }

    TEST(grating_test_1297){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 6, 6);
    }

    TEST(grating_test_1298){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 6, 7);
    }

    TEST(grating_test_1299){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 6, 7);
    }

    TEST(grating_test_1300){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 7, 1);
    }

    TEST(grating_test_1301){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 7, 1);
    }

    TEST(grating_test_1302){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 7, 2);
    }

    TEST(grating_test_1303){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 7, 2);
    }

    TEST(grating_test_1304){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 7, 3);
    }

    TEST(grating_test_1305){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 7, 3);
    }

    TEST(grating_test_1306){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 7, 5);
    }

    TEST(grating_test_1307){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 7, 5);
    }

    TEST(grating_test_1308){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 7, 6);
    }

    TEST(grating_test_1309){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 7, 6);
    }

    TEST(grating_test_1310){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 7, 7);
    }

    TEST(grating_test_1311){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 7, 7);
    }

    TEST(grating_test_1312){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 1, 1);
    }

    TEST(grating_test_1313){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 1, 1);
    }

    TEST(grating_test_1314){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 1, 2);
    }

    TEST(grating_test_1315){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 1, 2);
    }

    TEST(grating_test_1316){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 1, 3);
    }

    TEST(grating_test_1317){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 1, 3);
    }

    TEST(grating_test_1318){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 1, 5);
    }

    TEST(grating_test_1319){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 1, 5);
    }

    TEST(grating_test_1320){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 1, 6);
    }

    TEST(grating_test_1321){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 1, 6);
    }

    TEST(grating_test_1322){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 1, 7);
    }

    TEST(grating_test_1323){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 1, 7);
    }

    TEST(grating_test_1324){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 2, 1);
    }

    TEST(grating_test_1325){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 2, 1);
    }

    TEST(grating_test_1326){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 2, 2);
    }

    TEST(grating_test_1327){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 2, 2);
    }

    TEST(grating_test_1328){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 2, 3);
    }

    TEST(grating_test_1329){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 2, 3);
    }

    TEST(grating_test_1330){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 2, 5);
    }

    TEST(grating_test_1331){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 2, 5);
    }

    TEST(grating_test_1332){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 2, 6);
    }

    TEST(grating_test_1333){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 2, 6);
    }

    TEST(grating_test_1334){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 2, 7);
    }

    TEST(grating_test_1335){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 2, 7);
    }

    TEST(grating_test_1336){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 3, 1);
    }

    TEST(grating_test_1337){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 3, 1);
    }

    TEST(grating_test_1338){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 3, 2);
    }

    TEST(grating_test_1339){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 3, 2);
    }

    TEST(grating_test_1340){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 3, 3);
    }

    TEST(grating_test_1341){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 3, 3);
    }

    TEST(grating_test_1342){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 3, 5);
    }

    TEST(grating_test_1343){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 3, 5);
    }

    TEST(grating_test_1344){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 3, 6);
    }

    TEST(grating_test_1345){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 3, 6);
    }

    TEST(grating_test_1346){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 3, 7);
    }

    TEST(grating_test_1347){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 3, 7);
    }

    TEST(grating_test_1348){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 5, 1);
    }

    TEST(grating_test_1349){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 5, 1);
    }

    TEST(grating_test_1350){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 5, 2);
    }

    TEST(grating_test_1351){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 5, 2);
    }

    TEST(grating_test_1352){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 5, 3);
    }

    TEST(grating_test_1353){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 5, 3);
    }

    TEST(grating_test_1354){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 5, 5);
    }

    TEST(grating_test_1355){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 5, 5);
    }

    TEST(grating_test_1356){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 5, 6);
    }

    TEST(grating_test_1357){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 5, 6);
    }

    TEST(grating_test_1358){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 5, 7);
    }

    TEST(grating_test_1359){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 5, 7);
    }

    TEST(grating_test_1360){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 6, 1);
    }

    TEST(grating_test_1361){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 6, 1);
    }

    TEST(grating_test_1362){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 6, 2);
    }

    TEST(grating_test_1363){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 6, 2);
    }

    TEST(grating_test_1364){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 6, 3);
    }

    TEST(grating_test_1365){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 6, 3);
    }

    TEST(grating_test_1366){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 6, 5);
    }

    TEST(grating_test_1367){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 6, 5);
    }

    TEST(grating_test_1368){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 6, 6);
    }

    TEST(grating_test_1369){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 6, 6);
    }

    TEST(grating_test_1370){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 6, 7);
    }

    TEST(grating_test_1371){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 6, 7);
    }

    TEST(grating_test_1372){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 7, 1);
    }

    TEST(grating_test_1373){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 7, 1);
    }

    TEST(grating_test_1374){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 7, 2);
    }

    TEST(grating_test_1375){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 7, 2);
    }

    TEST(grating_test_1376){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 7, 3);
    }

    TEST(grating_test_1377){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 7, 3);
    }

    TEST(grating_test_1378){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 7, 5);
    }

    TEST(grating_test_1379){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 7, 5);
    }

    TEST(grating_test_1380){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 7, 6);
    }

    TEST(grating_test_1381){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 7, 6);
    }

    TEST(grating_test_1382){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 7, 7);
    }

    TEST(grating_test_1383){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 7, 7);
    }

    TEST(grating_test_1384){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 1, 1);
    }

    TEST(grating_test_1385){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 1, 1);
    }

    TEST(grating_test_1386){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 1, 2);
    }

    TEST(grating_test_1387){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 1, 2);
    }

    TEST(grating_test_1388){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 1, 3);
    }

    TEST(grating_test_1389){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 1, 3);
    }

    TEST(grating_test_1390){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 1, 5);
    }

    TEST(grating_test_1391){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 1, 5);
    }

    TEST(grating_test_1392){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 1, 6);
    }

    TEST(grating_test_1393){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 1, 6);
    }

    TEST(grating_test_1394){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 1, 7);
    }

    TEST(grating_test_1395){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 1, 7);
    }

    TEST(grating_test_1396){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 2, 1);
    }

    TEST(grating_test_1397){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 2, 1);
    }

    TEST(grating_test_1398){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 2, 2);
    }

    TEST(grating_test_1399){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 2, 2);
    }

    TEST(grating_test_1400){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 2, 3);
    }

    TEST(grating_test_1401){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 2, 3);
    }

    TEST(grating_test_1402){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 2, 5);
    }

    TEST(grating_test_1403){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 2, 5);
    }

    TEST(grating_test_1404){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 2, 6);
    }

    TEST(grating_test_1405){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 2, 6);
    }

    TEST(grating_test_1406){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 2, 7);
    }

    TEST(grating_test_1407){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 2, 7);
    }

    TEST(grating_test_1408){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 3, 1);
    }

    TEST(grating_test_1409){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 3, 1);
    }

    TEST(grating_test_1410){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 3, 2);
    }

    TEST(grating_test_1411){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 3, 2);
    }

    TEST(grating_test_1412){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 3, 3);
    }

    TEST(grating_test_1413){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 3, 3);
    }

    TEST(grating_test_1414){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 3, 5);
    }

    TEST(grating_test_1415){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 3, 5);
    }

    TEST(grating_test_1416){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 3, 6);
    }

    TEST(grating_test_1417){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 3, 6);
    }

    TEST(grating_test_1418){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 3, 7);
    }

    TEST(grating_test_1419){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 3, 7);
    }

    TEST(grating_test_1420){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 5, 1);
    }

    TEST(grating_test_1421){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 5, 1);
    }

    TEST(grating_test_1422){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 5, 2);
    }

    TEST(grating_test_1423){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 5, 2);
    }

    TEST(grating_test_1424){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 5, 3);
    }

    TEST(grating_test_1425){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 5, 3);
    }

    TEST(grating_test_1426){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 5, 5);
    }

    TEST(grating_test_1427){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 5, 5);
    }

    TEST(grating_test_1428){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 5, 6);
    }

    TEST(grating_test_1429){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 5, 6);
    }

    TEST(grating_test_1430){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 5, 7);
    }

    TEST(grating_test_1431){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 5, 7);
    }

    TEST(grating_test_1432){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 6, 1);
    }

    TEST(grating_test_1433){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 6, 1);
    }

    TEST(grating_test_1434){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 6, 2);
    }

    TEST(grating_test_1435){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 6, 2);
    }

    TEST(grating_test_1436){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 6, 3);
    }

    TEST(grating_test_1437){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 6, 3);
    }

    TEST(grating_test_1438){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 6, 5);
    }

    TEST(grating_test_1439){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 6, 5);
    }

    TEST(grating_test_1440){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 6, 6);
    }

    TEST(grating_test_1441){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 6, 6);
    }

    TEST(grating_test_1442){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 6, 7);
    }

    TEST(grating_test_1443){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 6, 7);
    }

    TEST(grating_test_1444){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 7, 1);
    }

    TEST(grating_test_1445){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 7, 1);
    }

    TEST(grating_test_1446){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 7, 2);
    }

    TEST(grating_test_1447){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 7, 2);
    }

    TEST(grating_test_1448){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 7, 3);
    }

    TEST(grating_test_1449){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 7, 3);
    }

    TEST(grating_test_1450){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 7, 5);
    }

    TEST(grating_test_1451){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 7, 5);
    }

    TEST(grating_test_1452){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 7, 6);
    }

    TEST(grating_test_1453){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 7, 6);
    }

    TEST(grating_test_1454){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 7, 7);
    }

    TEST(grating_test_1455){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 7, 7);
    }

    TEST(grating_test_1456){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 1, 1);
    }

    TEST(grating_test_1457){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 1, 1);
    }

    TEST(grating_test_1458){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 1, 2);
    }

    TEST(grating_test_1459){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 1, 2);
    }

    TEST(grating_test_1460){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 1, 3);
    }

    TEST(grating_test_1461){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 1, 3);
    }

    TEST(grating_test_1462){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 1, 5);
    }

    TEST(grating_test_1463){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 1, 5);
    }

    TEST(grating_test_1464){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 1, 6);
    }

    TEST(grating_test_1465){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 1, 6);
    }

    TEST(grating_test_1466){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 1, 7);
    }

    TEST(grating_test_1467){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 1, 7);
    }

    TEST(grating_test_1468){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 2, 1);
    }

    TEST(grating_test_1469){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 2, 1);
    }

    TEST(grating_test_1470){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 2, 2);
    }

    TEST(grating_test_1471){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 2, 2);
    }

    TEST(grating_test_1472){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 2, 3);
    }

    TEST(grating_test_1473){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 2, 3);
    }

    TEST(grating_test_1474){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 2, 5);
    }

    TEST(grating_test_1475){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 2, 5);
    }

    TEST(grating_test_1476){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 2, 6);
    }

    TEST(grating_test_1477){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 2, 6);
    }

    TEST(grating_test_1478){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 2, 7);
    }

    TEST(grating_test_1479){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 2, 7);
    }

    TEST(grating_test_1480){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 3, 1);
    }

    TEST(grating_test_1481){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 3, 1);
    }

    TEST(grating_test_1482){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 3, 2);
    }

    TEST(grating_test_1483){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 3, 2);
    }

    TEST(grating_test_1484){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 3, 3);
    }

    TEST(grating_test_1485){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 3, 3);
    }

    TEST(grating_test_1486){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 3, 5);
    }

    TEST(grating_test_1487){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 3, 5);
    }

    TEST(grating_test_1488){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 3, 6);
    }

    TEST(grating_test_1489){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 3, 6);
    }

    TEST(grating_test_1490){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 3, 7);
    }

    TEST(grating_test_1491){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 3, 7);
    }

    TEST(grating_test_1492){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 5, 1);
    }

    TEST(grating_test_1493){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 5, 1);
    }

    TEST(grating_test_1494){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 5, 2);
    }

    TEST(grating_test_1495){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 5, 2);
    }

    TEST(grating_test_1496){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 5, 3);
    }

    TEST(grating_test_1497){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 5, 3);
    }

    TEST(grating_test_1498){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 5, 5);
    }

    TEST(grating_test_1499){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 5, 5);
    }

    TEST(grating_test_1500){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 5, 6);
    }

    TEST(grating_test_1501){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 5, 6);
    }

    TEST(grating_test_1502){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 5, 7);
    }

    TEST(grating_test_1503){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 5, 7);
    }

    TEST(grating_test_1504){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 6, 1);
    }

    TEST(grating_test_1505){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 6, 1);
    }

    TEST(grating_test_1506){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 6, 2);
    }

    TEST(grating_test_1507){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 6, 2);
    }

    TEST(grating_test_1508){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 6, 3);
    }

    TEST(grating_test_1509){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 6, 3);
    }

    TEST(grating_test_1510){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 6, 5);
    }

    TEST(grating_test_1511){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 6, 5);
    }

    TEST(grating_test_1512){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 6, 6);
    }

    TEST(grating_test_1513){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 6, 6);
    }

    TEST(grating_test_1514){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 6, 7);
    }

    TEST(grating_test_1515){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 6, 7);
    }

    TEST(grating_test_1516){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 7, 1);
    }

    TEST(grating_test_1517){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 7, 1);
    }

    TEST(grating_test_1518){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 7, 2);
    }

    TEST(grating_test_1519){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 7, 2);
    }

    TEST(grating_test_1520){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 7, 3);
    }

    TEST(grating_test_1521){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 7, 3);
    }

    TEST(grating_test_1522){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 7, 5);
    }

    TEST(grating_test_1523){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 7, 5);
    }

    TEST(grating_test_1524){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 7, 6);
    }

    TEST(grating_test_1525){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 7, 6);
    }

    TEST(grating_test_1526){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 7, 7);
    }

    TEST(grating_test_1527){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 7, 7);
    }

    TEST(grating_test_1528){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 1, 1);
    }

    TEST(grating_test_1529){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 1, 1);
    }

    TEST(grating_test_1530){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 1, 2);
    }

    TEST(grating_test_1531){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 1, 2);
    }

    TEST(grating_test_1532){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 1, 3);
    }

    TEST(grating_test_1533){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 1, 3);
    }

    TEST(grating_test_1534){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 1, 5);
    }

    TEST(grating_test_1535){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 1, 5);
    }

    TEST(grating_test_1536){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 1, 6);
    }

    TEST(grating_test_1537){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 1, 6);
    }

    TEST(grating_test_1538){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 1, 7);
    }

    TEST(grating_test_1539){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 1, 7);
    }

    TEST(grating_test_1540){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 2, 1);
    }

    TEST(grating_test_1541){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 2, 1);
    }

    TEST(grating_test_1542){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 2, 2);
    }

    TEST(grating_test_1543){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 2, 2);
    }

    TEST(grating_test_1544){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 2, 3);
    }

    TEST(grating_test_1545){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 2, 3);
    }

    TEST(grating_test_1546){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 2, 5);
    }

    TEST(grating_test_1547){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 2, 5);
    }

    TEST(grating_test_1548){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 2, 6);
    }

    TEST(grating_test_1549){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 2, 6);
    }

    TEST(grating_test_1550){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 2, 7);
    }

    TEST(grating_test_1551){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 2, 7);
    }

    TEST(grating_test_1552){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 3, 1);
    }

    TEST(grating_test_1553){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 3, 1);
    }

    TEST(grating_test_1554){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 3, 2);
    }

    TEST(grating_test_1555){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 3, 2);
    }

    TEST(grating_test_1556){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 3, 3);
    }

    TEST(grating_test_1557){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 3, 3);
    }

    TEST(grating_test_1558){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 3, 5);
    }

    TEST(grating_test_1559){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 3, 5);
    }

    TEST(grating_test_1560){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 3, 6);
    }

    TEST(grating_test_1561){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 3, 6);
    }

    TEST(grating_test_1562){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 3, 7);
    }

    TEST(grating_test_1563){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 3, 7);
    }

    TEST(grating_test_1564){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 5, 1);
    }

    TEST(grating_test_1565){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 5, 1);
    }

    TEST(grating_test_1566){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 5, 2);
    }

    TEST(grating_test_1567){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 5, 2);
    }

    TEST(grating_test_1568){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 5, 3);
    }

    TEST(grating_test_1569){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 5, 3);
    }

    TEST(grating_test_1570){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 5, 5);
    }

    TEST(grating_test_1571){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 5, 5);
    }

    TEST(grating_test_1572){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 5, 6);
    }

    TEST(grating_test_1573){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 5, 6);
    }

    TEST(grating_test_1574){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 5, 7);
    }

    TEST(grating_test_1575){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 5, 7);
    }

    TEST(grating_test_1576){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 6, 1);
    }

    TEST(grating_test_1577){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 6, 1);
    }

    TEST(grating_test_1578){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 6, 2);
    }

    TEST(grating_test_1579){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 6, 2);
    }

    TEST(grating_test_1580){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 6, 3);
    }

    TEST(grating_test_1581){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 6, 3);
    }

    TEST(grating_test_1582){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 6, 5);
    }

    TEST(grating_test_1583){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 6, 5);
    }

    TEST(grating_test_1584){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 6, 6);
    }

    TEST(grating_test_1585){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 6, 6);
    }

    TEST(grating_test_1586){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 6, 7);
    }

    TEST(grating_test_1587){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 6, 7);
    }

    TEST(grating_test_1588){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 7, 1);
    }

    TEST(grating_test_1589){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 7, 1);
    }

    TEST(grating_test_1590){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 7, 2);
    }

    TEST(grating_test_1591){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 7, 2);
    }

    TEST(grating_test_1592){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 7, 3);
    }

    TEST(grating_test_1593){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 7, 3);
    }

    TEST(grating_test_1594){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 7, 5);
    }

    TEST(grating_test_1595){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 7, 5);
    }

    TEST(grating_test_1596){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 7, 6);
    }

    TEST(grating_test_1597){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 7, 6);
    }

    TEST(grating_test_1598){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 7, 7);
    }

    TEST(grating_test_1599){
         runIntegratorGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 7, 7);
    }


}


