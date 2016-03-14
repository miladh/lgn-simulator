/**********************************************************************
 *  Test: Full-field grating
 *
 *  Analytic source: closed-form experssion
 *
 * ********************************************************************/

#include <unittest++/UnitTest++.h>
#include <lgnSimulator.h>

using namespace lgnSimulator;

cube driftingGrating(double C, double wd, vec2 kd, vec r, vec t){

    cube S = zeros<cube>(r.n_elem, r.n_elem, t.n_elem);
    for(int l = 0; l < int(t.n_elem); l++){
        for(int i = 0; i < int(r.n_elem); i++){
            for(int j = 0; j < int(r.n_elem); j++){
                S(i,j,l) = C * cos(dot(kd, vec2{r[i], r[j]}) - wd * t[l]);
            }
        }
    }

    return S;
}

cx_cube driftingGratingFourierTransform(double C, double wd, vec2 kd, vec k, vec w){

    cx_cube S_ft = zeros<cx_cube>(k.n_elem, k.n_elem, w.n_elem);
    double dw = w[0] - w[1];
    double dk = k[1] - k[0];
    for(int l = 0; l < int(w.n_elem); l++){
        for(int i = 0; i < int(k.n_elem); i++){
            for(int j = 0; j < int(k.n_elem); j++){
                S_ft(i,j,l) = (Special::delta(vec2{k[i], k[j]}, kd)
                              * Special::delta(w[l], wd)
                            + Special::delta(vec2{k[i], k[j]}, -kd)
                               *Special::delta(-w[l], wd));
            }
        }
    }

    return S_ft * C * 4. * core::pi*core::pi*core::pi/dw/dk/dk;

}

void runStimulusFullFieldGratingTest(int ns, int nt, double dt, double ds,
                                     double C, int wdId, int kxId, int thetaId){


    Integrator integrator(nt, dt, ns, ds);
    vec r = integrator.spatialVec();
    vec k = integrator.spatialFreqVec();
    vec t = integrator.timeVec();
    vec w = integrator.temporalFreqVec();
    vec orientations = {0., 30., 45., -60., 90., -120., 180., -330.};

    double wd = w(wdId);
    double spatialFreq = k(kxId);
    double orientation = orientations(thetaId);


    FullFieldGrating grating(integrator, spatialFreq, orientation, wd, C);
    grating.computeFourierTransform();
    grating.computeSpatiotemporal();

    cube Sdg = driftingGrating(C, wd, grating.kVec(), r, t);
    cx_cube Sdg_ft = driftingGratingFourierTransform(C, wd ,grating.kVec(), k, w);

    cube diff_spatioTemporal = abs(grating.spatioTemporal() - Sdg);
    cx_cube diff_fourierTransform = grating.fourierTransform() - Sdg_ft;

    cube diff_fourierTransform_real = abs(real(diff_fourierTransform));
    cube diff_fourierTransform_imag = abs(imag(diff_fourierTransform));


    // Test
    CHECK_CLOSE(diff_spatioTemporal.max(), 0.0, 1e-9);
    CHECK_CLOSE(diff_fourierTransform_real.max(), 0.0, 1e-9);
    CHECK_CLOSE(diff_fourierTransform_imag.max(), 0.0, 1e-9);

}

SUITE(stimulus){


    TEST(FullFieldGrating_test_0){
         runStimulusFullFieldGratingTest(2, 2, 0.01, 0.01, -0.1, 0, 1, 1);
    }

    TEST(FullFieldGrating_test_1){
         runStimulusFullFieldGratingTest(2, 2, 0.01, 0.01, 2.2, 0, 1, 1);
    }

    TEST(FullFieldGrating_test_2){
         runStimulusFullFieldGratingTest(2, 2, 0.01, 0.01, -0.1, 0, 1, 3);
    }

    TEST(FullFieldGrating_test_3){
         runStimulusFullFieldGratingTest(2, 2, 0.01, 0.01, 2.2, 0, 1, 3);
    }

    TEST(FullFieldGrating_test_4){
         runStimulusFullFieldGratingTest(2, 2, 0.01, 0.01, -0.1, 0, 3, 1);
    }

    TEST(FullFieldGrating_test_5){
         runStimulusFullFieldGratingTest(2, 2, 0.01, 0.01, 2.2, 0, 3, 1);
    }

    TEST(FullFieldGrating_test_6){
         runStimulusFullFieldGratingTest(2, 2, 0.01, 0.01, -0.1, 0, 3, 3);
    }

    TEST(FullFieldGrating_test_7){
         runStimulusFullFieldGratingTest(2, 2, 0.01, 0.01, 2.2, 0, 3, 3);
    }

    TEST(FullFieldGrating_test_8){
         runStimulusFullFieldGratingTest(2, 2, 0.01, 0.01, -0.1, 1, 1, 1);
    }

    TEST(FullFieldGrating_test_9){
         runStimulusFullFieldGratingTest(2, 2, 0.01, 0.01, 2.2, 1, 1, 1);
    }

    TEST(FullFieldGrating_test_10){
         runStimulusFullFieldGratingTest(2, 2, 0.01, 0.01, -0.1, 1, 1, 3);
    }

    TEST(FullFieldGrating_test_11){
         runStimulusFullFieldGratingTest(2, 2, 0.01, 0.01, 2.2, 1, 1, 3);
    }

    TEST(FullFieldGrating_test_12){
         runStimulusFullFieldGratingTest(2, 2, 0.01, 0.01, -0.1, 1, 3, 1);
    }

    TEST(FullFieldGrating_test_13){
         runStimulusFullFieldGratingTest(2, 2, 0.01, 0.01, 2.2, 1, 3, 1);
    }

    TEST(FullFieldGrating_test_14){
         runStimulusFullFieldGratingTest(2, 2, 0.01, 0.01, -0.1, 1, 3, 3);
    }

    TEST(FullFieldGrating_test_15){
         runStimulusFullFieldGratingTest(2, 2, 0.01, 0.01, 2.2, 1, 3, 3);
    }

    TEST(FullFieldGrating_test_16){
         runStimulusFullFieldGratingTest(2, 2, 0.01, 0.01, -0.1, 3, 1, 1);
    }

    TEST(FullFieldGrating_test_17){
         runStimulusFullFieldGratingTest(2, 2, 0.01, 0.01, 2.2, 3, 1, 1);
    }

    TEST(FullFieldGrating_test_18){
         runStimulusFullFieldGratingTest(2, 2, 0.01, 0.01, -0.1, 3, 1, 3);
    }

    TEST(FullFieldGrating_test_19){
         runStimulusFullFieldGratingTest(2, 2, 0.01, 0.01, 2.2, 3, 1, 3);
    }

    TEST(FullFieldGrating_test_20){
         runStimulusFullFieldGratingTest(2, 2, 0.01, 0.01, -0.1, 3, 3, 1);
    }

    TEST(FullFieldGrating_test_21){
         runStimulusFullFieldGratingTest(2, 2, 0.01, 0.01, 2.2, 3, 3, 1);
    }

    TEST(FullFieldGrating_test_22){
         runStimulusFullFieldGratingTest(2, 2, 0.01, 0.01, -0.1, 3, 3, 3);
    }

    TEST(FullFieldGrating_test_23){
         runStimulusFullFieldGratingTest(2, 2, 0.01, 0.01, 2.2, 3, 3, 3);
    }

    TEST(FullFieldGrating_test_24){
         runStimulusFullFieldGratingTest(2, 2, 0.5, 0.01, -0.1, 0, 1, 1);
    }

    TEST(FullFieldGrating_test_25){
         runStimulusFullFieldGratingTest(2, 2, 0.5, 0.01, 2.2, 0, 1, 1);
    }

    TEST(FullFieldGrating_test_26){
         runStimulusFullFieldGratingTest(2, 2, 0.5, 0.01, -0.1, 0, 1, 3);
    }

    TEST(FullFieldGrating_test_27){
         runStimulusFullFieldGratingTest(2, 2, 0.5, 0.01, 2.2, 0, 1, 3);
    }

    TEST(FullFieldGrating_test_28){
         runStimulusFullFieldGratingTest(2, 2, 0.5, 0.01, -0.1, 0, 3, 1);
    }

    TEST(FullFieldGrating_test_29){
         runStimulusFullFieldGratingTest(2, 2, 0.5, 0.01, 2.2, 0, 3, 1);
    }

    TEST(FullFieldGrating_test_30){
         runStimulusFullFieldGratingTest(2, 2, 0.5, 0.01, -0.1, 0, 3, 3);
    }

    TEST(FullFieldGrating_test_31){
         runStimulusFullFieldGratingTest(2, 2, 0.5, 0.01, 2.2, 0, 3, 3);
    }

    TEST(FullFieldGrating_test_32){
         runStimulusFullFieldGratingTest(2, 2, 0.5, 0.01, -0.1, 1, 1, 1);
    }

    TEST(FullFieldGrating_test_33){
         runStimulusFullFieldGratingTest(2, 2, 0.5, 0.01, 2.2, 1, 1, 1);
    }

    TEST(FullFieldGrating_test_34){
         runStimulusFullFieldGratingTest(2, 2, 0.5, 0.01, -0.1, 1, 1, 3);
    }

    TEST(FullFieldGrating_test_35){
         runStimulusFullFieldGratingTest(2, 2, 0.5, 0.01, 2.2, 1, 1, 3);
    }

    TEST(FullFieldGrating_test_36){
         runStimulusFullFieldGratingTest(2, 2, 0.5, 0.01, -0.1, 1, 3, 1);
    }

    TEST(FullFieldGrating_test_37){
         runStimulusFullFieldGratingTest(2, 2, 0.5, 0.01, 2.2, 1, 3, 1);
    }

    TEST(FullFieldGrating_test_38){
         runStimulusFullFieldGratingTest(2, 2, 0.5, 0.01, -0.1, 1, 3, 3);
    }

    TEST(FullFieldGrating_test_39){
         runStimulusFullFieldGratingTest(2, 2, 0.5, 0.01, 2.2, 1, 3, 3);
    }

    TEST(FullFieldGrating_test_40){
         runStimulusFullFieldGratingTest(2, 2, 0.5, 0.01, -0.1, 3, 1, 1);
    }

    TEST(FullFieldGrating_test_41){
         runStimulusFullFieldGratingTest(2, 2, 0.5, 0.01, 2.2, 3, 1, 1);
    }

    TEST(FullFieldGrating_test_42){
         runStimulusFullFieldGratingTest(2, 2, 0.5, 0.01, -0.1, 3, 1, 3);
    }

    TEST(FullFieldGrating_test_43){
         runStimulusFullFieldGratingTest(2, 2, 0.5, 0.01, 2.2, 3, 1, 3);
    }

    TEST(FullFieldGrating_test_44){
         runStimulusFullFieldGratingTest(2, 2, 0.5, 0.01, -0.1, 3, 3, 1);
    }

    TEST(FullFieldGrating_test_45){
         runStimulusFullFieldGratingTest(2, 2, 0.5, 0.01, 2.2, 3, 3, 1);
    }

    TEST(FullFieldGrating_test_46){
         runStimulusFullFieldGratingTest(2, 2, 0.5, 0.01, -0.1, 3, 3, 3);
    }

    TEST(FullFieldGrating_test_47){
         runStimulusFullFieldGratingTest(2, 2, 0.5, 0.01, 2.2, 3, 3, 3);
    }

    TEST(FullFieldGrating_test_48){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, -0.1, 0, 1, 1);
    }

    TEST(FullFieldGrating_test_49){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, 2.2, 0, 1, 1);
    }

    TEST(FullFieldGrating_test_50){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, -0.1, 0, 1, 3);
    }

    TEST(FullFieldGrating_test_51){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, 2.2, 0, 1, 3);
    }

    TEST(FullFieldGrating_test_52){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, -0.1, 0, 3, 1);
    }

    TEST(FullFieldGrating_test_53){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, 2.2, 0, 3, 1);
    }

    TEST(FullFieldGrating_test_54){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, -0.1, 0, 3, 3);
    }

    TEST(FullFieldGrating_test_55){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, 2.2, 0, 3, 3);
    }

    TEST(FullFieldGrating_test_56){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, -0.1, 1, 1, 1);
    }

    TEST(FullFieldGrating_test_57){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, 2.2, 1, 1, 1);
    }

    TEST(FullFieldGrating_test_58){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, -0.1, 1, 1, 3);
    }

    TEST(FullFieldGrating_test_59){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, 2.2, 1, 1, 3);
    }

    TEST(FullFieldGrating_test_60){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, -0.1, 1, 3, 1);
    }

    TEST(FullFieldGrating_test_61){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, 2.2, 1, 3, 1);
    }

    TEST(FullFieldGrating_test_62){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, -0.1, 1, 3, 3);
    }

    TEST(FullFieldGrating_test_63){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, 2.2, 1, 3, 3);
    }

    TEST(FullFieldGrating_test_64){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, -0.1, 2, 1, 1);
    }

    TEST(FullFieldGrating_test_65){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, 2.2, 2, 1, 1);
    }

    TEST(FullFieldGrating_test_66){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, -0.1, 2, 1, 3);
    }

    TEST(FullFieldGrating_test_67){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, 2.2, 2, 1, 3);
    }

    TEST(FullFieldGrating_test_68){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, -0.1, 2, 3, 1);
    }

    TEST(FullFieldGrating_test_69){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, 2.2, 2, 3, 1);
    }

    TEST(FullFieldGrating_test_70){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, -0.1, 2, 3, 3);
    }

    TEST(FullFieldGrating_test_71){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, 2.2, 2, 3, 3);
    }

    TEST(FullFieldGrating_test_72){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, -0.1, 3, 1, 1);
    }

    TEST(FullFieldGrating_test_73){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, 2.2, 3, 1, 1);
    }

    TEST(FullFieldGrating_test_74){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, -0.1, 3, 1, 3);
    }

    TEST(FullFieldGrating_test_75){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, 2.2, 3, 1, 3);
    }

    TEST(FullFieldGrating_test_76){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, -0.1, 3, 3, 1);
    }

    TEST(FullFieldGrating_test_77){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, 2.2, 3, 3, 1);
    }

    TEST(FullFieldGrating_test_78){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, -0.1, 3, 3, 3);
    }

    TEST(FullFieldGrating_test_79){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, 2.2, 3, 3, 3);
    }

    TEST(FullFieldGrating_test_80){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, -0.1, 5, 1, 1);
    }

    TEST(FullFieldGrating_test_81){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, 2.2, 5, 1, 1);
    }

    TEST(FullFieldGrating_test_82){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, -0.1, 5, 1, 3);
    }

    TEST(FullFieldGrating_test_83){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, 2.2, 5, 1, 3);
    }

    TEST(FullFieldGrating_test_84){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, -0.1, 5, 3, 1);
    }

    TEST(FullFieldGrating_test_85){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, 2.2, 5, 3, 1);
    }

    TEST(FullFieldGrating_test_86){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, -0.1, 5, 3, 3);
    }

    TEST(FullFieldGrating_test_87){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, 2.2, 5, 3, 3);
    }

    TEST(FullFieldGrating_test_88){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, -0.1, 6, 1, 1);
    }

    TEST(FullFieldGrating_test_89){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, 2.2, 6, 1, 1);
    }

    TEST(FullFieldGrating_test_90){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, -0.1, 6, 1, 3);
    }

    TEST(FullFieldGrating_test_91){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, 2.2, 6, 1, 3);
    }

    TEST(FullFieldGrating_test_92){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, -0.1, 6, 3, 1);
    }

    TEST(FullFieldGrating_test_93){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, 2.2, 6, 3, 1);
    }

    TEST(FullFieldGrating_test_94){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, -0.1, 6, 3, 3);
    }

    TEST(FullFieldGrating_test_95){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, 2.2, 6, 3, 3);
    }

    TEST(FullFieldGrating_test_96){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, -0.1, 7, 1, 1);
    }

    TEST(FullFieldGrating_test_97){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, 2.2, 7, 1, 1);
    }

    TEST(FullFieldGrating_test_98){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, -0.1, 7, 1, 3);
    }

    TEST(FullFieldGrating_test_99){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, 2.2, 7, 1, 3);
    }

    TEST(FullFieldGrating_test_100){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, -0.1, 7, 3, 1);
    }

    TEST(FullFieldGrating_test_101){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, 2.2, 7, 3, 1);
    }

    TEST(FullFieldGrating_test_102){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, -0.1, 7, 3, 3);
    }

    TEST(FullFieldGrating_test_103){
         runStimulusFullFieldGratingTest(2, 3, 0.01, 0.01, 2.2, 7, 3, 3);
    }

    TEST(FullFieldGrating_test_104){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, -0.1, 0, 1, 1);
    }

    TEST(FullFieldGrating_test_105){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, 2.2, 0, 1, 1);
    }

    TEST(FullFieldGrating_test_106){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, -0.1, 0, 1, 3);
    }

    TEST(FullFieldGrating_test_107){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, 2.2, 0, 1, 3);
    }

    TEST(FullFieldGrating_test_108){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, -0.1, 0, 3, 1);
    }

    TEST(FullFieldGrating_test_109){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, 2.2, 0, 3, 1);
    }

    TEST(FullFieldGrating_test_110){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, -0.1, 0, 3, 3);
    }

    TEST(FullFieldGrating_test_111){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, 2.2, 0, 3, 3);
    }

    TEST(FullFieldGrating_test_112){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, -0.1, 1, 1, 1);
    }

    TEST(FullFieldGrating_test_113){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, 2.2, 1, 1, 1);
    }

    TEST(FullFieldGrating_test_114){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, -0.1, 1, 1, 3);
    }

    TEST(FullFieldGrating_test_115){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, 2.2, 1, 1, 3);
    }

    TEST(FullFieldGrating_test_116){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, -0.1, 1, 3, 1);
    }

    TEST(FullFieldGrating_test_117){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, 2.2, 1, 3, 1);
    }

    TEST(FullFieldGrating_test_118){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, -0.1, 1, 3, 3);
    }

    TEST(FullFieldGrating_test_119){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, 2.2, 1, 3, 3);
    }

    TEST(FullFieldGrating_test_120){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, -0.1, 2, 1, 1);
    }

    TEST(FullFieldGrating_test_121){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, 2.2, 2, 1, 1);
    }

    TEST(FullFieldGrating_test_122){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, -0.1, 2, 1, 3);
    }

    TEST(FullFieldGrating_test_123){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, 2.2, 2, 1, 3);
    }

    TEST(FullFieldGrating_test_124){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, -0.1, 2, 3, 1);
    }

    TEST(FullFieldGrating_test_125){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, 2.2, 2, 3, 1);
    }

    TEST(FullFieldGrating_test_126){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, -0.1, 2, 3, 3);
    }

    TEST(FullFieldGrating_test_127){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, 2.2, 2, 3, 3);
    }

    TEST(FullFieldGrating_test_128){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, -0.1, 3, 1, 1);
    }

    TEST(FullFieldGrating_test_129){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, 2.2, 3, 1, 1);
    }

    TEST(FullFieldGrating_test_130){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, -0.1, 3, 1, 3);
    }

    TEST(FullFieldGrating_test_131){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, 2.2, 3, 1, 3);
    }

    TEST(FullFieldGrating_test_132){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, -0.1, 3, 3, 1);
    }

    TEST(FullFieldGrating_test_133){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, 2.2, 3, 3, 1);
    }

    TEST(FullFieldGrating_test_134){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, -0.1, 3, 3, 3);
    }

    TEST(FullFieldGrating_test_135){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, 2.2, 3, 3, 3);
    }

    TEST(FullFieldGrating_test_136){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, -0.1, 5, 1, 1);
    }

    TEST(FullFieldGrating_test_137){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, 2.2, 5, 1, 1);
    }

    TEST(FullFieldGrating_test_138){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, -0.1, 5, 1, 3);
    }

    TEST(FullFieldGrating_test_139){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, 2.2, 5, 1, 3);
    }

    TEST(FullFieldGrating_test_140){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, -0.1, 5, 3, 1);
    }

    TEST(FullFieldGrating_test_141){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, 2.2, 5, 3, 1);
    }

    TEST(FullFieldGrating_test_142){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, -0.1, 5, 3, 3);
    }

    TEST(FullFieldGrating_test_143){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, 2.2, 5, 3, 3);
    }

    TEST(FullFieldGrating_test_144){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, -0.1, 6, 1, 1);
    }

    TEST(FullFieldGrating_test_145){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, 2.2, 6, 1, 1);
    }

    TEST(FullFieldGrating_test_146){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, -0.1, 6, 1, 3);
    }

    TEST(FullFieldGrating_test_147){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, 2.2, 6, 1, 3);
    }

    TEST(FullFieldGrating_test_148){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, -0.1, 6, 3, 1);
    }

    TEST(FullFieldGrating_test_149){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, 2.2, 6, 3, 1);
    }

    TEST(FullFieldGrating_test_150){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, -0.1, 6, 3, 3);
    }

    TEST(FullFieldGrating_test_151){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, 2.2, 6, 3, 3);
    }

    TEST(FullFieldGrating_test_152){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, -0.1, 7, 1, 1);
    }

    TEST(FullFieldGrating_test_153){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, 2.2, 7, 1, 1);
    }

    TEST(FullFieldGrating_test_154){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, -0.1, 7, 1, 3);
    }

    TEST(FullFieldGrating_test_155){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, 2.2, 7, 1, 3);
    }

    TEST(FullFieldGrating_test_156){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, -0.1, 7, 3, 1);
    }

    TEST(FullFieldGrating_test_157){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, 2.2, 7, 3, 1);
    }

    TEST(FullFieldGrating_test_158){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, -0.1, 7, 3, 3);
    }

    TEST(FullFieldGrating_test_159){
         runStimulusFullFieldGratingTest(2, 3, 0.5, 0.01, 2.2, 7, 3, 3);
    }

    TEST(FullFieldGrating_test_160){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 1, 1);
    }

    TEST(FullFieldGrating_test_161){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 1, 1);
    }

    TEST(FullFieldGrating_test_162){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 1, 2);
    }

    TEST(FullFieldGrating_test_163){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 1, 2);
    }

    TEST(FullFieldGrating_test_164){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 1, 3);
    }

    TEST(FullFieldGrating_test_165){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 1, 3);
    }

    TEST(FullFieldGrating_test_166){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 1, 5);
    }

    TEST(FullFieldGrating_test_167){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 1, 5);
    }

    TEST(FullFieldGrating_test_168){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 1, 6);
    }

    TEST(FullFieldGrating_test_169){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 1, 6);
    }

    TEST(FullFieldGrating_test_170){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 1, 7);
    }

    TEST(FullFieldGrating_test_171){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 1, 7);
    }

    TEST(FullFieldGrating_test_172){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 2, 1);
    }

    TEST(FullFieldGrating_test_173){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 2, 1);
    }

    TEST(FullFieldGrating_test_174){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 2, 2);
    }

    TEST(FullFieldGrating_test_175){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 2, 2);
    }

    TEST(FullFieldGrating_test_176){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 2, 3);
    }

    TEST(FullFieldGrating_test_177){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 2, 3);
    }

    TEST(FullFieldGrating_test_178){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 2, 5);
    }

    TEST(FullFieldGrating_test_179){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 2, 5);
    }

    TEST(FullFieldGrating_test_180){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 2, 6);
    }

    TEST(FullFieldGrating_test_181){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 2, 6);
    }

    TEST(FullFieldGrating_test_182){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 2, 7);
    }

    TEST(FullFieldGrating_test_183){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 2, 7);
    }

    TEST(FullFieldGrating_test_184){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 3, 1);
    }

    TEST(FullFieldGrating_test_185){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 3, 1);
    }

    TEST(FullFieldGrating_test_186){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 3, 2);
    }

    TEST(FullFieldGrating_test_187){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 3, 2);
    }

    TEST(FullFieldGrating_test_188){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 3, 3);
    }

    TEST(FullFieldGrating_test_189){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 3, 3);
    }

    TEST(FullFieldGrating_test_190){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 3, 5);
    }

    TEST(FullFieldGrating_test_191){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 3, 5);
    }

    TEST(FullFieldGrating_test_192){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 3, 6);
    }

    TEST(FullFieldGrating_test_193){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 3, 6);
    }

    TEST(FullFieldGrating_test_194){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 3, 7);
    }

    TEST(FullFieldGrating_test_195){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 3, 7);
    }

    TEST(FullFieldGrating_test_196){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 5, 1);
    }

    TEST(FullFieldGrating_test_197){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 5, 1);
    }

    TEST(FullFieldGrating_test_198){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 5, 2);
    }

    TEST(FullFieldGrating_test_199){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 5, 2);
    }

    TEST(FullFieldGrating_test_200){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 5, 3);
    }

    TEST(FullFieldGrating_test_201){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 5, 3);
    }

    TEST(FullFieldGrating_test_202){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 5, 5);
    }

    TEST(FullFieldGrating_test_203){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 5, 5);
    }

    TEST(FullFieldGrating_test_204){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 5, 6);
    }

    TEST(FullFieldGrating_test_205){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 5, 6);
    }

    TEST(FullFieldGrating_test_206){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 5, 7);
    }

    TEST(FullFieldGrating_test_207){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 5, 7);
    }

    TEST(FullFieldGrating_test_208){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 6, 1);
    }

    TEST(FullFieldGrating_test_209){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 6, 1);
    }

    TEST(FullFieldGrating_test_210){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 6, 2);
    }

    TEST(FullFieldGrating_test_211){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 6, 2);
    }

    TEST(FullFieldGrating_test_212){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 6, 3);
    }

    TEST(FullFieldGrating_test_213){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 6, 3);
    }

    TEST(FullFieldGrating_test_214){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 6, 5);
    }

    TEST(FullFieldGrating_test_215){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 6, 5);
    }

    TEST(FullFieldGrating_test_216){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 6, 6);
    }

    TEST(FullFieldGrating_test_217){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 6, 6);
    }

    TEST(FullFieldGrating_test_218){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 6, 7);
    }

    TEST(FullFieldGrating_test_219){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 6, 7);
    }

    TEST(FullFieldGrating_test_220){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 7, 1);
    }

    TEST(FullFieldGrating_test_221){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 7, 1);
    }

    TEST(FullFieldGrating_test_222){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 7, 2);
    }

    TEST(FullFieldGrating_test_223){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 7, 2);
    }

    TEST(FullFieldGrating_test_224){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 7, 3);
    }

    TEST(FullFieldGrating_test_225){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 7, 3);
    }

    TEST(FullFieldGrating_test_226){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 7, 5);
    }

    TEST(FullFieldGrating_test_227){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 7, 5);
    }

    TEST(FullFieldGrating_test_228){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 7, 6);
    }

    TEST(FullFieldGrating_test_229){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 7, 6);
    }

    TEST(FullFieldGrating_test_230){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 0, 7, 7);
    }

    TEST(FullFieldGrating_test_231){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 0, 7, 7);
    }

    TEST(FullFieldGrating_test_232){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 1, 1);
    }

    TEST(FullFieldGrating_test_233){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 1, 1);
    }

    TEST(FullFieldGrating_test_234){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 1, 2);
    }

    TEST(FullFieldGrating_test_235){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 1, 2);
    }

    TEST(FullFieldGrating_test_236){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 1, 3);
    }

    TEST(FullFieldGrating_test_237){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 1, 3);
    }

    TEST(FullFieldGrating_test_238){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 1, 5);
    }

    TEST(FullFieldGrating_test_239){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 1, 5);
    }

    TEST(FullFieldGrating_test_240){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 1, 6);
    }

    TEST(FullFieldGrating_test_241){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 1, 6);
    }

    TEST(FullFieldGrating_test_242){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 1, 7);
    }

    TEST(FullFieldGrating_test_243){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 1, 7);
    }

    TEST(FullFieldGrating_test_244){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 2, 1);
    }

    TEST(FullFieldGrating_test_245){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 2, 1);
    }

    TEST(FullFieldGrating_test_246){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 2, 2);
    }

    TEST(FullFieldGrating_test_247){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 2, 2);
    }

    TEST(FullFieldGrating_test_248){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 2, 3);
    }

    TEST(FullFieldGrating_test_249){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 2, 3);
    }

    TEST(FullFieldGrating_test_250){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 2, 5);
    }

    TEST(FullFieldGrating_test_251){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 2, 5);
    }

    TEST(FullFieldGrating_test_252){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 2, 6);
    }

    TEST(FullFieldGrating_test_253){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 2, 6);
    }

    TEST(FullFieldGrating_test_254){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 2, 7);
    }

    TEST(FullFieldGrating_test_255){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 2, 7);
    }

    TEST(FullFieldGrating_test_256){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 3, 1);
    }

    TEST(FullFieldGrating_test_257){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 3, 1);
    }

    TEST(FullFieldGrating_test_258){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 3, 2);
    }

    TEST(FullFieldGrating_test_259){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 3, 2);
    }

    TEST(FullFieldGrating_test_260){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 3, 3);
    }

    TEST(FullFieldGrating_test_261){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 3, 3);
    }

    TEST(FullFieldGrating_test_262){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 3, 5);
    }

    TEST(FullFieldGrating_test_263){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 3, 5);
    }

    TEST(FullFieldGrating_test_264){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 3, 6);
    }

    TEST(FullFieldGrating_test_265){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 3, 6);
    }

    TEST(FullFieldGrating_test_266){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 3, 7);
    }

    TEST(FullFieldGrating_test_267){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 3, 7);
    }

    TEST(FullFieldGrating_test_268){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 5, 1);
    }

    TEST(FullFieldGrating_test_269){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 5, 1);
    }

    TEST(FullFieldGrating_test_270){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 5, 2);
    }

    TEST(FullFieldGrating_test_271){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 5, 2);
    }

    TEST(FullFieldGrating_test_272){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 5, 3);
    }

    TEST(FullFieldGrating_test_273){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 5, 3);
    }

    TEST(FullFieldGrating_test_274){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 5, 5);
    }

    TEST(FullFieldGrating_test_275){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 5, 5);
    }

    TEST(FullFieldGrating_test_276){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 5, 6);
    }

    TEST(FullFieldGrating_test_277){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 5, 6);
    }

    TEST(FullFieldGrating_test_278){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 5, 7);
    }

    TEST(FullFieldGrating_test_279){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 5, 7);
    }

    TEST(FullFieldGrating_test_280){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 6, 1);
    }

    TEST(FullFieldGrating_test_281){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 6, 1);
    }

    TEST(FullFieldGrating_test_282){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 6, 2);
    }

    TEST(FullFieldGrating_test_283){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 6, 2);
    }

    TEST(FullFieldGrating_test_284){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 6, 3);
    }

    TEST(FullFieldGrating_test_285){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 6, 3);
    }

    TEST(FullFieldGrating_test_286){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 6, 5);
    }

    TEST(FullFieldGrating_test_287){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 6, 5);
    }

    TEST(FullFieldGrating_test_288){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 6, 6);
    }

    TEST(FullFieldGrating_test_289){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 6, 6);
    }

    TEST(FullFieldGrating_test_290){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 6, 7);
    }

    TEST(FullFieldGrating_test_291){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 6, 7);
    }

    TEST(FullFieldGrating_test_292){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 7, 1);
    }

    TEST(FullFieldGrating_test_293){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 7, 1);
    }

    TEST(FullFieldGrating_test_294){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 7, 2);
    }

    TEST(FullFieldGrating_test_295){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 7, 2);
    }

    TEST(FullFieldGrating_test_296){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 7, 3);
    }

    TEST(FullFieldGrating_test_297){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 7, 3);
    }

    TEST(FullFieldGrating_test_298){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 7, 5);
    }

    TEST(FullFieldGrating_test_299){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 7, 5);
    }

    TEST(FullFieldGrating_test_300){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 7, 6);
    }

    TEST(FullFieldGrating_test_301){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 7, 6);
    }

    TEST(FullFieldGrating_test_302){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 1, 7, 7);
    }

    TEST(FullFieldGrating_test_303){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 1, 7, 7);
    }

    TEST(FullFieldGrating_test_304){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 1, 1);
    }

    TEST(FullFieldGrating_test_305){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 1, 1);
    }

    TEST(FullFieldGrating_test_306){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 1, 2);
    }

    TEST(FullFieldGrating_test_307){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 1, 2);
    }

    TEST(FullFieldGrating_test_308){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 1, 3);
    }

    TEST(FullFieldGrating_test_309){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 1, 3);
    }

    TEST(FullFieldGrating_test_310){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 1, 5);
    }

    TEST(FullFieldGrating_test_311){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 1, 5);
    }

    TEST(FullFieldGrating_test_312){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 1, 6);
    }

    TEST(FullFieldGrating_test_313){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 1, 6);
    }

    TEST(FullFieldGrating_test_314){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 1, 7);
    }

    TEST(FullFieldGrating_test_315){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 1, 7);
    }

    TEST(FullFieldGrating_test_316){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 2, 1);
    }

    TEST(FullFieldGrating_test_317){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 2, 1);
    }

    TEST(FullFieldGrating_test_318){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 2, 2);
    }

    TEST(FullFieldGrating_test_319){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 2, 2);
    }

    TEST(FullFieldGrating_test_320){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 2, 3);
    }

    TEST(FullFieldGrating_test_321){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 2, 3);
    }

    TEST(FullFieldGrating_test_322){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 2, 5);
    }

    TEST(FullFieldGrating_test_323){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 2, 5);
    }

    TEST(FullFieldGrating_test_324){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 2, 6);
    }

    TEST(FullFieldGrating_test_325){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 2, 6);
    }

    TEST(FullFieldGrating_test_326){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 2, 7);
    }

    TEST(FullFieldGrating_test_327){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 2, 7);
    }

    TEST(FullFieldGrating_test_328){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 3, 1);
    }

    TEST(FullFieldGrating_test_329){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 3, 1);
    }

    TEST(FullFieldGrating_test_330){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 3, 2);
    }

    TEST(FullFieldGrating_test_331){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 3, 2);
    }

    TEST(FullFieldGrating_test_332){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 3, 3);
    }

    TEST(FullFieldGrating_test_333){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 3, 3);
    }

    TEST(FullFieldGrating_test_334){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 3, 5);
    }

    TEST(FullFieldGrating_test_335){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 3, 5);
    }

    TEST(FullFieldGrating_test_336){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 3, 6);
    }

    TEST(FullFieldGrating_test_337){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 3, 6);
    }

    TEST(FullFieldGrating_test_338){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 3, 7);
    }

    TEST(FullFieldGrating_test_339){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 3, 7);
    }

    TEST(FullFieldGrating_test_340){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 5, 1);
    }

    TEST(FullFieldGrating_test_341){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 5, 1);
    }

    TEST(FullFieldGrating_test_342){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 5, 2);
    }

    TEST(FullFieldGrating_test_343){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 5, 2);
    }

    TEST(FullFieldGrating_test_344){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 5, 3);
    }

    TEST(FullFieldGrating_test_345){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 5, 3);
    }

    TEST(FullFieldGrating_test_346){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 5, 5);
    }

    TEST(FullFieldGrating_test_347){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 5, 5);
    }

    TEST(FullFieldGrating_test_348){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 5, 6);
    }

    TEST(FullFieldGrating_test_349){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 5, 6);
    }

    TEST(FullFieldGrating_test_350){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 5, 7);
    }

    TEST(FullFieldGrating_test_351){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 5, 7);
    }

    TEST(FullFieldGrating_test_352){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 6, 1);
    }

    TEST(FullFieldGrating_test_353){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 6, 1);
    }

    TEST(FullFieldGrating_test_354){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 6, 2);
    }

    TEST(FullFieldGrating_test_355){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 6, 2);
    }

    TEST(FullFieldGrating_test_356){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 6, 3);
    }

    TEST(FullFieldGrating_test_357){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 6, 3);
    }

    TEST(FullFieldGrating_test_358){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 6, 5);
    }

    TEST(FullFieldGrating_test_359){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 6, 5);
    }

    TEST(FullFieldGrating_test_360){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 6, 6);
    }

    TEST(FullFieldGrating_test_361){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 6, 6);
    }

    TEST(FullFieldGrating_test_362){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 6, 7);
    }

    TEST(FullFieldGrating_test_363){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 6, 7);
    }

    TEST(FullFieldGrating_test_364){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 7, 1);
    }

    TEST(FullFieldGrating_test_365){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 7, 1);
    }

    TEST(FullFieldGrating_test_366){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 7, 2);
    }

    TEST(FullFieldGrating_test_367){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 7, 2);
    }

    TEST(FullFieldGrating_test_368){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 7, 3);
    }

    TEST(FullFieldGrating_test_369){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 7, 3);
    }

    TEST(FullFieldGrating_test_370){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 7, 5);
    }

    TEST(FullFieldGrating_test_371){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 7, 5);
    }

    TEST(FullFieldGrating_test_372){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 7, 6);
    }

    TEST(FullFieldGrating_test_373){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 7, 6);
    }

    TEST(FullFieldGrating_test_374){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, -0.1, 3, 7, 7);
    }

    TEST(FullFieldGrating_test_375){
         runStimulusFullFieldGratingTest(3, 2, 0.01, 0.01, 2.2, 3, 7, 7);
    }

    TEST(FullFieldGrating_test_376){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 1, 1);
    }

    TEST(FullFieldGrating_test_377){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 1, 1);
    }

    TEST(FullFieldGrating_test_378){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 1, 2);
    }

    TEST(FullFieldGrating_test_379){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 1, 2);
    }

    TEST(FullFieldGrating_test_380){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 1, 3);
    }

    TEST(FullFieldGrating_test_381){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 1, 3);
    }

    TEST(FullFieldGrating_test_382){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 1, 5);
    }

    TEST(FullFieldGrating_test_383){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 1, 5);
    }

    TEST(FullFieldGrating_test_384){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 1, 6);
    }

    TEST(FullFieldGrating_test_385){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 1, 6);
    }

    TEST(FullFieldGrating_test_386){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 1, 7);
    }

    TEST(FullFieldGrating_test_387){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 1, 7);
    }

    TEST(FullFieldGrating_test_388){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 2, 1);
    }

    TEST(FullFieldGrating_test_389){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 2, 1);
    }

    TEST(FullFieldGrating_test_390){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 2, 2);
    }

    TEST(FullFieldGrating_test_391){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 2, 2);
    }

    TEST(FullFieldGrating_test_392){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 2, 3);
    }

    TEST(FullFieldGrating_test_393){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 2, 3);
    }

    TEST(FullFieldGrating_test_394){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 2, 5);
    }

    TEST(FullFieldGrating_test_395){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 2, 5);
    }

    TEST(FullFieldGrating_test_396){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 2, 6);
    }

    TEST(FullFieldGrating_test_397){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 2, 6);
    }

    TEST(FullFieldGrating_test_398){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 2, 7);
    }

    TEST(FullFieldGrating_test_399){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 2, 7);
    }

    TEST(FullFieldGrating_test_400){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 3, 1);
    }

    TEST(FullFieldGrating_test_401){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 3, 1);
    }

    TEST(FullFieldGrating_test_402){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 3, 2);
    }

    TEST(FullFieldGrating_test_403){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 3, 2);
    }

    TEST(FullFieldGrating_test_404){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 3, 3);
    }

    TEST(FullFieldGrating_test_405){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 3, 3);
    }

    TEST(FullFieldGrating_test_406){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 3, 5);
    }

    TEST(FullFieldGrating_test_407){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 3, 5);
    }

    TEST(FullFieldGrating_test_408){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 3, 6);
    }

    TEST(FullFieldGrating_test_409){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 3, 6);
    }

    TEST(FullFieldGrating_test_410){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 3, 7);
    }

    TEST(FullFieldGrating_test_411){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 3, 7);
    }

    TEST(FullFieldGrating_test_412){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 5, 1);
    }

    TEST(FullFieldGrating_test_413){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 5, 1);
    }

    TEST(FullFieldGrating_test_414){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 5, 2);
    }

    TEST(FullFieldGrating_test_415){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 5, 2);
    }

    TEST(FullFieldGrating_test_416){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 5, 3);
    }

    TEST(FullFieldGrating_test_417){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 5, 3);
    }

    TEST(FullFieldGrating_test_418){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 5, 5);
    }

    TEST(FullFieldGrating_test_419){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 5, 5);
    }

    TEST(FullFieldGrating_test_420){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 5, 6);
    }

    TEST(FullFieldGrating_test_421){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 5, 6);
    }

    TEST(FullFieldGrating_test_422){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 5, 7);
    }

    TEST(FullFieldGrating_test_423){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 5, 7);
    }

    TEST(FullFieldGrating_test_424){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 6, 1);
    }

    TEST(FullFieldGrating_test_425){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 6, 1);
    }

    TEST(FullFieldGrating_test_426){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 6, 2);
    }

    TEST(FullFieldGrating_test_427){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 6, 2);
    }

    TEST(FullFieldGrating_test_428){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 6, 3);
    }

    TEST(FullFieldGrating_test_429){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 6, 3);
    }

    TEST(FullFieldGrating_test_430){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 6, 5);
    }

    TEST(FullFieldGrating_test_431){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 6, 5);
    }

    TEST(FullFieldGrating_test_432){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 6, 6);
    }

    TEST(FullFieldGrating_test_433){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 6, 6);
    }

    TEST(FullFieldGrating_test_434){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 6, 7);
    }

    TEST(FullFieldGrating_test_435){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 6, 7);
    }

    TEST(FullFieldGrating_test_436){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 7, 1);
    }

    TEST(FullFieldGrating_test_437){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 7, 1);
    }

    TEST(FullFieldGrating_test_438){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 7, 2);
    }

    TEST(FullFieldGrating_test_439){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 7, 2);
    }

    TEST(FullFieldGrating_test_440){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 7, 3);
    }

    TEST(FullFieldGrating_test_441){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 7, 3);
    }

    TEST(FullFieldGrating_test_442){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 7, 5);
    }

    TEST(FullFieldGrating_test_443){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 7, 5);
    }

    TEST(FullFieldGrating_test_444){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 7, 6);
    }

    TEST(FullFieldGrating_test_445){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 7, 6);
    }

    TEST(FullFieldGrating_test_446){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 0, 7, 7);
    }

    TEST(FullFieldGrating_test_447){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 0, 7, 7);
    }

    TEST(FullFieldGrating_test_448){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 1, 1);
    }

    TEST(FullFieldGrating_test_449){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 1, 1);
    }

    TEST(FullFieldGrating_test_450){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 1, 2);
    }

    TEST(FullFieldGrating_test_451){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 1, 2);
    }

    TEST(FullFieldGrating_test_452){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 1, 3);
    }

    TEST(FullFieldGrating_test_453){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 1, 3);
    }

    TEST(FullFieldGrating_test_454){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 1, 5);
    }

    TEST(FullFieldGrating_test_455){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 1, 5);
    }

    TEST(FullFieldGrating_test_456){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 1, 6);
    }

    TEST(FullFieldGrating_test_457){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 1, 6);
    }

    TEST(FullFieldGrating_test_458){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 1, 7);
    }

    TEST(FullFieldGrating_test_459){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 1, 7);
    }

    TEST(FullFieldGrating_test_460){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 2, 1);
    }

    TEST(FullFieldGrating_test_461){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 2, 1);
    }

    TEST(FullFieldGrating_test_462){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 2, 2);
    }

    TEST(FullFieldGrating_test_463){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 2, 2);
    }

    TEST(FullFieldGrating_test_464){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 2, 3);
    }

    TEST(FullFieldGrating_test_465){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 2, 3);
    }

    TEST(FullFieldGrating_test_466){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 2, 5);
    }

    TEST(FullFieldGrating_test_467){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 2, 5);
    }

    TEST(FullFieldGrating_test_468){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 2, 6);
    }

    TEST(FullFieldGrating_test_469){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 2, 6);
    }

    TEST(FullFieldGrating_test_470){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 2, 7);
    }

    TEST(FullFieldGrating_test_471){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 2, 7);
    }

    TEST(FullFieldGrating_test_472){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 3, 1);
    }

    TEST(FullFieldGrating_test_473){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 3, 1);
    }

    TEST(FullFieldGrating_test_474){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 3, 2);
    }

    TEST(FullFieldGrating_test_475){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 3, 2);
    }

    TEST(FullFieldGrating_test_476){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 3, 3);
    }

    TEST(FullFieldGrating_test_477){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 3, 3);
    }

    TEST(FullFieldGrating_test_478){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 3, 5);
    }

    TEST(FullFieldGrating_test_479){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 3, 5);
    }

    TEST(FullFieldGrating_test_480){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 3, 6);
    }

    TEST(FullFieldGrating_test_481){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 3, 6);
    }

    TEST(FullFieldGrating_test_482){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 3, 7);
    }

    TEST(FullFieldGrating_test_483){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 3, 7);
    }

    TEST(FullFieldGrating_test_484){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 5, 1);
    }

    TEST(FullFieldGrating_test_485){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 5, 1);
    }

    TEST(FullFieldGrating_test_486){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 5, 2);
    }

    TEST(FullFieldGrating_test_487){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 5, 2);
    }

    TEST(FullFieldGrating_test_488){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 5, 3);
    }

    TEST(FullFieldGrating_test_489){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 5, 3);
    }

    TEST(FullFieldGrating_test_490){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 5, 5);
    }

    TEST(FullFieldGrating_test_491){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 5, 5);
    }

    TEST(FullFieldGrating_test_492){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 5, 6);
    }

    TEST(FullFieldGrating_test_493){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 5, 6);
    }

    TEST(FullFieldGrating_test_494){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 5, 7);
    }

    TEST(FullFieldGrating_test_495){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 5, 7);
    }

    TEST(FullFieldGrating_test_496){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 6, 1);
    }

    TEST(FullFieldGrating_test_497){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 6, 1);
    }

    TEST(FullFieldGrating_test_498){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 6, 2);
    }

    TEST(FullFieldGrating_test_499){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 6, 2);
    }

    TEST(FullFieldGrating_test_500){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 6, 3);
    }

    TEST(FullFieldGrating_test_501){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 6, 3);
    }

    TEST(FullFieldGrating_test_502){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 6, 5);
    }

    TEST(FullFieldGrating_test_503){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 6, 5);
    }

    TEST(FullFieldGrating_test_504){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 6, 6);
    }

    TEST(FullFieldGrating_test_505){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 6, 6);
    }

    TEST(FullFieldGrating_test_506){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 6, 7);
    }

    TEST(FullFieldGrating_test_507){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 6, 7);
    }

    TEST(FullFieldGrating_test_508){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 7, 1);
    }

    TEST(FullFieldGrating_test_509){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 7, 1);
    }

    TEST(FullFieldGrating_test_510){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 7, 2);
    }

    TEST(FullFieldGrating_test_511){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 7, 2);
    }

    TEST(FullFieldGrating_test_512){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 7, 3);
    }

    TEST(FullFieldGrating_test_513){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 7, 3);
    }

    TEST(FullFieldGrating_test_514){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 7, 5);
    }

    TEST(FullFieldGrating_test_515){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 7, 5);
    }

    TEST(FullFieldGrating_test_516){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 7, 6);
    }

    TEST(FullFieldGrating_test_517){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 7, 6);
    }

    TEST(FullFieldGrating_test_518){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 1, 7, 7);
    }

    TEST(FullFieldGrating_test_519){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 1, 7, 7);
    }

    TEST(FullFieldGrating_test_520){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 1, 1);
    }

    TEST(FullFieldGrating_test_521){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 1, 1);
    }

    TEST(FullFieldGrating_test_522){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 1, 2);
    }

    TEST(FullFieldGrating_test_523){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 1, 2);
    }

    TEST(FullFieldGrating_test_524){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 1, 3);
    }

    TEST(FullFieldGrating_test_525){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 1, 3);
    }

    TEST(FullFieldGrating_test_526){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 1, 5);
    }

    TEST(FullFieldGrating_test_527){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 1, 5);
    }

    TEST(FullFieldGrating_test_528){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 1, 6);
    }

    TEST(FullFieldGrating_test_529){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 1, 6);
    }

    TEST(FullFieldGrating_test_530){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 1, 7);
    }

    TEST(FullFieldGrating_test_531){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 1, 7);
    }

    TEST(FullFieldGrating_test_532){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 2, 1);
    }

    TEST(FullFieldGrating_test_533){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 2, 1);
    }

    TEST(FullFieldGrating_test_534){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 2, 2);
    }

    TEST(FullFieldGrating_test_535){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 2, 2);
    }

    TEST(FullFieldGrating_test_536){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 2, 3);
    }

    TEST(FullFieldGrating_test_537){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 2, 3);
    }

    TEST(FullFieldGrating_test_538){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 2, 5);
    }

    TEST(FullFieldGrating_test_539){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 2, 5);
    }

    TEST(FullFieldGrating_test_540){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 2, 6);
    }

    TEST(FullFieldGrating_test_541){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 2, 6);
    }

    TEST(FullFieldGrating_test_542){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 2, 7);
    }

    TEST(FullFieldGrating_test_543){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 2, 7);
    }

    TEST(FullFieldGrating_test_544){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 3, 1);
    }

    TEST(FullFieldGrating_test_545){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 3, 1);
    }

    TEST(FullFieldGrating_test_546){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 3, 2);
    }

    TEST(FullFieldGrating_test_547){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 3, 2);
    }

    TEST(FullFieldGrating_test_548){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 3, 3);
    }

    TEST(FullFieldGrating_test_549){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 3, 3);
    }

    TEST(FullFieldGrating_test_550){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 3, 5);
    }

    TEST(FullFieldGrating_test_551){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 3, 5);
    }

    TEST(FullFieldGrating_test_552){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 3, 6);
    }

    TEST(FullFieldGrating_test_553){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 3, 6);
    }

    TEST(FullFieldGrating_test_554){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 3, 7);
    }

    TEST(FullFieldGrating_test_555){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 3, 7);
    }

    TEST(FullFieldGrating_test_556){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 5, 1);
    }

    TEST(FullFieldGrating_test_557){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 5, 1);
    }

    TEST(FullFieldGrating_test_558){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 5, 2);
    }

    TEST(FullFieldGrating_test_559){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 5, 2);
    }

    TEST(FullFieldGrating_test_560){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 5, 3);
    }

    TEST(FullFieldGrating_test_561){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 5, 3);
    }

    TEST(FullFieldGrating_test_562){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 5, 5);
    }

    TEST(FullFieldGrating_test_563){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 5, 5);
    }

    TEST(FullFieldGrating_test_564){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 5, 6);
    }

    TEST(FullFieldGrating_test_565){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 5, 6);
    }

    TEST(FullFieldGrating_test_566){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 5, 7);
    }

    TEST(FullFieldGrating_test_567){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 5, 7);
    }

    TEST(FullFieldGrating_test_568){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 6, 1);
    }

    TEST(FullFieldGrating_test_569){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 6, 1);
    }

    TEST(FullFieldGrating_test_570){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 6, 2);
    }

    TEST(FullFieldGrating_test_571){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 6, 2);
    }

    TEST(FullFieldGrating_test_572){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 6, 3);
    }

    TEST(FullFieldGrating_test_573){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 6, 3);
    }

    TEST(FullFieldGrating_test_574){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 6, 5);
    }

    TEST(FullFieldGrating_test_575){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 6, 5);
    }

    TEST(FullFieldGrating_test_576){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 6, 6);
    }

    TEST(FullFieldGrating_test_577){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 6, 6);
    }

    TEST(FullFieldGrating_test_578){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 6, 7);
    }

    TEST(FullFieldGrating_test_579){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 6, 7);
    }

    TEST(FullFieldGrating_test_580){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 7, 1);
    }

    TEST(FullFieldGrating_test_581){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 7, 1);
    }

    TEST(FullFieldGrating_test_582){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 7, 2);
    }

    TEST(FullFieldGrating_test_583){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 7, 2);
    }

    TEST(FullFieldGrating_test_584){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 7, 3);
    }

    TEST(FullFieldGrating_test_585){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 7, 3);
    }

    TEST(FullFieldGrating_test_586){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 7, 5);
    }

    TEST(FullFieldGrating_test_587){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 7, 5);
    }

    TEST(FullFieldGrating_test_588){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 7, 6);
    }

    TEST(FullFieldGrating_test_589){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 7, 6);
    }

    TEST(FullFieldGrating_test_590){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, -0.1, 3, 7, 7);
    }

    TEST(FullFieldGrating_test_591){
         runStimulusFullFieldGratingTest(3, 2, 0.5, 0.01, 2.2, 3, 7, 7);
    }

    TEST(FullFieldGrating_test_592){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 1, 1);
    }

    TEST(FullFieldGrating_test_593){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 1, 1);
    }

    TEST(FullFieldGrating_test_594){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 1, 2);
    }

    TEST(FullFieldGrating_test_595){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 1, 2);
    }

    TEST(FullFieldGrating_test_596){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 1, 3);
    }

    TEST(FullFieldGrating_test_597){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 1, 3);
    }

    TEST(FullFieldGrating_test_598){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 1, 5);
    }

    TEST(FullFieldGrating_test_599){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 1, 5);
    }

    TEST(FullFieldGrating_test_600){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 1, 6);
    }

    TEST(FullFieldGrating_test_601){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 1, 6);
    }

    TEST(FullFieldGrating_test_602){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 1, 7);
    }

    TEST(FullFieldGrating_test_603){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 1, 7);
    }

    TEST(FullFieldGrating_test_604){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 2, 1);
    }

    TEST(FullFieldGrating_test_605){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 2, 1);
    }

    TEST(FullFieldGrating_test_606){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 2, 2);
    }

    TEST(FullFieldGrating_test_607){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 2, 2);
    }

    TEST(FullFieldGrating_test_608){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 2, 3);
    }

    TEST(FullFieldGrating_test_609){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 2, 3);
    }

    TEST(FullFieldGrating_test_610){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 2, 5);
    }

    TEST(FullFieldGrating_test_611){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 2, 5);
    }

    TEST(FullFieldGrating_test_612){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 2, 6);
    }

    TEST(FullFieldGrating_test_613){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 2, 6);
    }

    TEST(FullFieldGrating_test_614){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 2, 7);
    }

    TEST(FullFieldGrating_test_615){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 2, 7);
    }

    TEST(FullFieldGrating_test_616){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 3, 1);
    }

    TEST(FullFieldGrating_test_617){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 3, 1);
    }

    TEST(FullFieldGrating_test_618){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 3, 2);
    }

    TEST(FullFieldGrating_test_619){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 3, 2);
    }

    TEST(FullFieldGrating_test_620){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 3, 3);
    }

    TEST(FullFieldGrating_test_621){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 3, 3);
    }

    TEST(FullFieldGrating_test_622){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 3, 5);
    }

    TEST(FullFieldGrating_test_623){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 3, 5);
    }

    TEST(FullFieldGrating_test_624){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 3, 6);
    }

    TEST(FullFieldGrating_test_625){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 3, 6);
    }

    TEST(FullFieldGrating_test_626){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 3, 7);
    }

    TEST(FullFieldGrating_test_627){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 3, 7);
    }

    TEST(FullFieldGrating_test_628){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 5, 1);
    }

    TEST(FullFieldGrating_test_629){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 5, 1);
    }

    TEST(FullFieldGrating_test_630){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 5, 2);
    }

    TEST(FullFieldGrating_test_631){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 5, 2);
    }

    TEST(FullFieldGrating_test_632){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 5, 3);
    }

    TEST(FullFieldGrating_test_633){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 5, 3);
    }

    TEST(FullFieldGrating_test_634){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 5, 5);
    }

    TEST(FullFieldGrating_test_635){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 5, 5);
    }

    TEST(FullFieldGrating_test_636){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 5, 6);
    }

    TEST(FullFieldGrating_test_637){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 5, 6);
    }

    TEST(FullFieldGrating_test_638){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 5, 7);
    }

    TEST(FullFieldGrating_test_639){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 5, 7);
    }

    TEST(FullFieldGrating_test_640){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 6, 1);
    }

    TEST(FullFieldGrating_test_641){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 6, 1);
    }

    TEST(FullFieldGrating_test_642){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 6, 2);
    }

    TEST(FullFieldGrating_test_643){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 6, 2);
    }

    TEST(FullFieldGrating_test_644){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 6, 3);
    }

    TEST(FullFieldGrating_test_645){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 6, 3);
    }

    TEST(FullFieldGrating_test_646){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 6, 5);
    }

    TEST(FullFieldGrating_test_647){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 6, 5);
    }

    TEST(FullFieldGrating_test_648){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 6, 6);
    }

    TEST(FullFieldGrating_test_649){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 6, 6);
    }

    TEST(FullFieldGrating_test_650){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 6, 7);
    }

    TEST(FullFieldGrating_test_651){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 6, 7);
    }

    TEST(FullFieldGrating_test_652){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 7, 1);
    }

    TEST(FullFieldGrating_test_653){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 7, 1);
    }

    TEST(FullFieldGrating_test_654){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 7, 2);
    }

    TEST(FullFieldGrating_test_655){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 7, 2);
    }

    TEST(FullFieldGrating_test_656){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 7, 3);
    }

    TEST(FullFieldGrating_test_657){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 7, 3);
    }

    TEST(FullFieldGrating_test_658){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 7, 5);
    }

    TEST(FullFieldGrating_test_659){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 7, 5);
    }

    TEST(FullFieldGrating_test_660){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 7, 6);
    }

    TEST(FullFieldGrating_test_661){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 7, 6);
    }

    TEST(FullFieldGrating_test_662){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 0, 7, 7);
    }

    TEST(FullFieldGrating_test_663){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 0, 7, 7);
    }

    TEST(FullFieldGrating_test_664){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 1, 1);
    }

    TEST(FullFieldGrating_test_665){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 1, 1);
    }

    TEST(FullFieldGrating_test_666){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 1, 2);
    }

    TEST(FullFieldGrating_test_667){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 1, 2);
    }

    TEST(FullFieldGrating_test_668){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 1, 3);
    }

    TEST(FullFieldGrating_test_669){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 1, 3);
    }

    TEST(FullFieldGrating_test_670){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 1, 5);
    }

    TEST(FullFieldGrating_test_671){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 1, 5);
    }

    TEST(FullFieldGrating_test_672){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 1, 6);
    }

    TEST(FullFieldGrating_test_673){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 1, 6);
    }

    TEST(FullFieldGrating_test_674){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 1, 7);
    }

    TEST(FullFieldGrating_test_675){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 1, 7);
    }

    TEST(FullFieldGrating_test_676){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 2, 1);
    }

    TEST(FullFieldGrating_test_677){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 2, 1);
    }

    TEST(FullFieldGrating_test_678){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 2, 2);
    }

    TEST(FullFieldGrating_test_679){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 2, 2);
    }

    TEST(FullFieldGrating_test_680){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 2, 3);
    }

    TEST(FullFieldGrating_test_681){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 2, 3);
    }

    TEST(FullFieldGrating_test_682){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 2, 5);
    }

    TEST(FullFieldGrating_test_683){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 2, 5);
    }

    TEST(FullFieldGrating_test_684){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 2, 6);
    }

    TEST(FullFieldGrating_test_685){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 2, 6);
    }

    TEST(FullFieldGrating_test_686){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 2, 7);
    }

    TEST(FullFieldGrating_test_687){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 2, 7);
    }

    TEST(FullFieldGrating_test_688){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 3, 1);
    }

    TEST(FullFieldGrating_test_689){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 3, 1);
    }

    TEST(FullFieldGrating_test_690){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 3, 2);
    }

    TEST(FullFieldGrating_test_691){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 3, 2);
    }

    TEST(FullFieldGrating_test_692){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 3, 3);
    }

    TEST(FullFieldGrating_test_693){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 3, 3);
    }

    TEST(FullFieldGrating_test_694){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 3, 5);
    }

    TEST(FullFieldGrating_test_695){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 3, 5);
    }

    TEST(FullFieldGrating_test_696){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 3, 6);
    }

    TEST(FullFieldGrating_test_697){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 3, 6);
    }

    TEST(FullFieldGrating_test_698){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 3, 7);
    }

    TEST(FullFieldGrating_test_699){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 3, 7);
    }

    TEST(FullFieldGrating_test_700){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 5, 1);
    }

    TEST(FullFieldGrating_test_701){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 5, 1);
    }

    TEST(FullFieldGrating_test_702){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 5, 2);
    }

    TEST(FullFieldGrating_test_703){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 5, 2);
    }

    TEST(FullFieldGrating_test_704){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 5, 3);
    }

    TEST(FullFieldGrating_test_705){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 5, 3);
    }

    TEST(FullFieldGrating_test_706){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 5, 5);
    }

    TEST(FullFieldGrating_test_707){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 5, 5);
    }

    TEST(FullFieldGrating_test_708){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 5, 6);
    }

    TEST(FullFieldGrating_test_709){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 5, 6);
    }

    TEST(FullFieldGrating_test_710){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 5, 7);
    }

    TEST(FullFieldGrating_test_711){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 5, 7);
    }

    TEST(FullFieldGrating_test_712){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 6, 1);
    }

    TEST(FullFieldGrating_test_713){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 6, 1);
    }

    TEST(FullFieldGrating_test_714){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 6, 2);
    }

    TEST(FullFieldGrating_test_715){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 6, 2);
    }

    TEST(FullFieldGrating_test_716){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 6, 3);
    }

    TEST(FullFieldGrating_test_717){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 6, 3);
    }

    TEST(FullFieldGrating_test_718){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 6, 5);
    }

    TEST(FullFieldGrating_test_719){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 6, 5);
    }

    TEST(FullFieldGrating_test_720){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 6, 6);
    }

    TEST(FullFieldGrating_test_721){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 6, 6);
    }

    TEST(FullFieldGrating_test_722){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 6, 7);
    }

    TEST(FullFieldGrating_test_723){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 6, 7);
    }

    TEST(FullFieldGrating_test_724){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 7, 1);
    }

    TEST(FullFieldGrating_test_725){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 7, 1);
    }

    TEST(FullFieldGrating_test_726){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 7, 2);
    }

    TEST(FullFieldGrating_test_727){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 7, 2);
    }

    TEST(FullFieldGrating_test_728){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 7, 3);
    }

    TEST(FullFieldGrating_test_729){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 7, 3);
    }

    TEST(FullFieldGrating_test_730){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 7, 5);
    }

    TEST(FullFieldGrating_test_731){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 7, 5);
    }

    TEST(FullFieldGrating_test_732){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 7, 6);
    }

    TEST(FullFieldGrating_test_733){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 7, 6);
    }

    TEST(FullFieldGrating_test_734){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 1, 7, 7);
    }

    TEST(FullFieldGrating_test_735){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 1, 7, 7);
    }

    TEST(FullFieldGrating_test_736){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 1, 1);
    }

    TEST(FullFieldGrating_test_737){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 1, 1);
    }

    TEST(FullFieldGrating_test_738){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 1, 2);
    }

    TEST(FullFieldGrating_test_739){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 1, 2);
    }

    TEST(FullFieldGrating_test_740){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 1, 3);
    }

    TEST(FullFieldGrating_test_741){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 1, 3);
    }

    TEST(FullFieldGrating_test_742){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 1, 5);
    }

    TEST(FullFieldGrating_test_743){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 1, 5);
    }

    TEST(FullFieldGrating_test_744){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 1, 6);
    }

    TEST(FullFieldGrating_test_745){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 1, 6);
    }

    TEST(FullFieldGrating_test_746){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 1, 7);
    }

    TEST(FullFieldGrating_test_747){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 1, 7);
    }

    TEST(FullFieldGrating_test_748){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 2, 1);
    }

    TEST(FullFieldGrating_test_749){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 2, 1);
    }

    TEST(FullFieldGrating_test_750){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 2, 2);
    }

    TEST(FullFieldGrating_test_751){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 2, 2);
    }

    TEST(FullFieldGrating_test_752){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 2, 3);
    }

    TEST(FullFieldGrating_test_753){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 2, 3);
    }

    TEST(FullFieldGrating_test_754){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 2, 5);
    }

    TEST(FullFieldGrating_test_755){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 2, 5);
    }

    TEST(FullFieldGrating_test_756){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 2, 6);
    }

    TEST(FullFieldGrating_test_757){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 2, 6);
    }

    TEST(FullFieldGrating_test_758){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 2, 7);
    }

    TEST(FullFieldGrating_test_759){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 2, 7);
    }

    TEST(FullFieldGrating_test_760){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 3, 1);
    }

    TEST(FullFieldGrating_test_761){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 3, 1);
    }

    TEST(FullFieldGrating_test_762){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 3, 2);
    }

    TEST(FullFieldGrating_test_763){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 3, 2);
    }

    TEST(FullFieldGrating_test_764){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 3, 3);
    }

    TEST(FullFieldGrating_test_765){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 3, 3);
    }

    TEST(FullFieldGrating_test_766){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 3, 5);
    }

    TEST(FullFieldGrating_test_767){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 3, 5);
    }

    TEST(FullFieldGrating_test_768){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 3, 6);
    }

    TEST(FullFieldGrating_test_769){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 3, 6);
    }

    TEST(FullFieldGrating_test_770){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 3, 7);
    }

    TEST(FullFieldGrating_test_771){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 3, 7);
    }

    TEST(FullFieldGrating_test_772){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 5, 1);
    }

    TEST(FullFieldGrating_test_773){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 5, 1);
    }

    TEST(FullFieldGrating_test_774){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 5, 2);
    }

    TEST(FullFieldGrating_test_775){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 5, 2);
    }

    TEST(FullFieldGrating_test_776){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 5, 3);
    }

    TEST(FullFieldGrating_test_777){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 5, 3);
    }

    TEST(FullFieldGrating_test_778){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 5, 5);
    }

    TEST(FullFieldGrating_test_779){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 5, 5);
    }

    TEST(FullFieldGrating_test_780){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 5, 6);
    }

    TEST(FullFieldGrating_test_781){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 5, 6);
    }

    TEST(FullFieldGrating_test_782){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 5, 7);
    }

    TEST(FullFieldGrating_test_783){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 5, 7);
    }

    TEST(FullFieldGrating_test_784){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 6, 1);
    }

    TEST(FullFieldGrating_test_785){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 6, 1);
    }

    TEST(FullFieldGrating_test_786){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 6, 2);
    }

    TEST(FullFieldGrating_test_787){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 6, 2);
    }

    TEST(FullFieldGrating_test_788){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 6, 3);
    }

    TEST(FullFieldGrating_test_789){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 6, 3);
    }

    TEST(FullFieldGrating_test_790){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 6, 5);
    }

    TEST(FullFieldGrating_test_791){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 6, 5);
    }

    TEST(FullFieldGrating_test_792){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 6, 6);
    }

    TEST(FullFieldGrating_test_793){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 6, 6);
    }

    TEST(FullFieldGrating_test_794){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 6, 7);
    }

    TEST(FullFieldGrating_test_795){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 6, 7);
    }

    TEST(FullFieldGrating_test_796){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 7, 1);
    }

    TEST(FullFieldGrating_test_797){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 7, 1);
    }

    TEST(FullFieldGrating_test_798){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 7, 2);
    }

    TEST(FullFieldGrating_test_799){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 7, 2);
    }

    TEST(FullFieldGrating_test_800){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 7, 3);
    }

    TEST(FullFieldGrating_test_801){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 7, 3);
    }

    TEST(FullFieldGrating_test_802){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 7, 5);
    }

    TEST(FullFieldGrating_test_803){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 7, 5);
    }

    TEST(FullFieldGrating_test_804){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 7, 6);
    }

    TEST(FullFieldGrating_test_805){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 7, 6);
    }

    TEST(FullFieldGrating_test_806){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 2, 7, 7);
    }

    TEST(FullFieldGrating_test_807){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 2, 7, 7);
    }

    TEST(FullFieldGrating_test_808){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 1, 1);
    }

    TEST(FullFieldGrating_test_809){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 1, 1);
    }

    TEST(FullFieldGrating_test_810){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 1, 2);
    }

    TEST(FullFieldGrating_test_811){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 1, 2);
    }

    TEST(FullFieldGrating_test_812){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 1, 3);
    }

    TEST(FullFieldGrating_test_813){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 1, 3);
    }

    TEST(FullFieldGrating_test_814){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 1, 5);
    }

    TEST(FullFieldGrating_test_815){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 1, 5);
    }

    TEST(FullFieldGrating_test_816){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 1, 6);
    }

    TEST(FullFieldGrating_test_817){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 1, 6);
    }

    TEST(FullFieldGrating_test_818){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 1, 7);
    }

    TEST(FullFieldGrating_test_819){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 1, 7);
    }

    TEST(FullFieldGrating_test_820){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 2, 1);
    }

    TEST(FullFieldGrating_test_821){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 2, 1);
    }

    TEST(FullFieldGrating_test_822){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 2, 2);
    }

    TEST(FullFieldGrating_test_823){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 2, 2);
    }

    TEST(FullFieldGrating_test_824){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 2, 3);
    }

    TEST(FullFieldGrating_test_825){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 2, 3);
    }

    TEST(FullFieldGrating_test_826){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 2, 5);
    }

    TEST(FullFieldGrating_test_827){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 2, 5);
    }

    TEST(FullFieldGrating_test_828){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 2, 6);
    }

    TEST(FullFieldGrating_test_829){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 2, 6);
    }

    TEST(FullFieldGrating_test_830){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 2, 7);
    }

    TEST(FullFieldGrating_test_831){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 2, 7);
    }

    TEST(FullFieldGrating_test_832){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 3, 1);
    }

    TEST(FullFieldGrating_test_833){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 3, 1);
    }

    TEST(FullFieldGrating_test_834){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 3, 2);
    }

    TEST(FullFieldGrating_test_835){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 3, 2);
    }

    TEST(FullFieldGrating_test_836){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 3, 3);
    }

    TEST(FullFieldGrating_test_837){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 3, 3);
    }

    TEST(FullFieldGrating_test_838){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 3, 5);
    }

    TEST(FullFieldGrating_test_839){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 3, 5);
    }

    TEST(FullFieldGrating_test_840){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 3, 6);
    }

    TEST(FullFieldGrating_test_841){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 3, 6);
    }

    TEST(FullFieldGrating_test_842){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 3, 7);
    }

    TEST(FullFieldGrating_test_843){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 3, 7);
    }

    TEST(FullFieldGrating_test_844){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 5, 1);
    }

    TEST(FullFieldGrating_test_845){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 5, 1);
    }

    TEST(FullFieldGrating_test_846){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 5, 2);
    }

    TEST(FullFieldGrating_test_847){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 5, 2);
    }

    TEST(FullFieldGrating_test_848){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 5, 3);
    }

    TEST(FullFieldGrating_test_849){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 5, 3);
    }

    TEST(FullFieldGrating_test_850){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 5, 5);
    }

    TEST(FullFieldGrating_test_851){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 5, 5);
    }

    TEST(FullFieldGrating_test_852){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 5, 6);
    }

    TEST(FullFieldGrating_test_853){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 5, 6);
    }

    TEST(FullFieldGrating_test_854){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 5, 7);
    }

    TEST(FullFieldGrating_test_855){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 5, 7);
    }

    TEST(FullFieldGrating_test_856){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 6, 1);
    }

    TEST(FullFieldGrating_test_857){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 6, 1);
    }

    TEST(FullFieldGrating_test_858){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 6, 2);
    }

    TEST(FullFieldGrating_test_859){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 6, 2);
    }

    TEST(FullFieldGrating_test_860){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 6, 3);
    }

    TEST(FullFieldGrating_test_861){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 6, 3);
    }

    TEST(FullFieldGrating_test_862){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 6, 5);
    }

    TEST(FullFieldGrating_test_863){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 6, 5);
    }

    TEST(FullFieldGrating_test_864){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 6, 6);
    }

    TEST(FullFieldGrating_test_865){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 6, 6);
    }

    TEST(FullFieldGrating_test_866){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 6, 7);
    }

    TEST(FullFieldGrating_test_867){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 6, 7);
    }

    TEST(FullFieldGrating_test_868){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 7, 1);
    }

    TEST(FullFieldGrating_test_869){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 7, 1);
    }

    TEST(FullFieldGrating_test_870){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 7, 2);
    }

    TEST(FullFieldGrating_test_871){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 7, 2);
    }

    TEST(FullFieldGrating_test_872){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 7, 3);
    }

    TEST(FullFieldGrating_test_873){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 7, 3);
    }

    TEST(FullFieldGrating_test_874){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 7, 5);
    }

    TEST(FullFieldGrating_test_875){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 7, 5);
    }

    TEST(FullFieldGrating_test_876){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 7, 6);
    }

    TEST(FullFieldGrating_test_877){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 7, 6);
    }

    TEST(FullFieldGrating_test_878){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 3, 7, 7);
    }

    TEST(FullFieldGrating_test_879){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 3, 7, 7);
    }

    TEST(FullFieldGrating_test_880){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 1, 1);
    }

    TEST(FullFieldGrating_test_881){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 1, 1);
    }

    TEST(FullFieldGrating_test_882){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 1, 2);
    }

    TEST(FullFieldGrating_test_883){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 1, 2);
    }

    TEST(FullFieldGrating_test_884){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 1, 3);
    }

    TEST(FullFieldGrating_test_885){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 1, 3);
    }

    TEST(FullFieldGrating_test_886){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 1, 5);
    }

    TEST(FullFieldGrating_test_887){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 1, 5);
    }

    TEST(FullFieldGrating_test_888){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 1, 6);
    }

    TEST(FullFieldGrating_test_889){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 1, 6);
    }

    TEST(FullFieldGrating_test_890){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 1, 7);
    }

    TEST(FullFieldGrating_test_891){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 1, 7);
    }

    TEST(FullFieldGrating_test_892){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 2, 1);
    }

    TEST(FullFieldGrating_test_893){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 2, 1);
    }

    TEST(FullFieldGrating_test_894){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 2, 2);
    }

    TEST(FullFieldGrating_test_895){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 2, 2);
    }

    TEST(FullFieldGrating_test_896){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 2, 3);
    }

    TEST(FullFieldGrating_test_897){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 2, 3);
    }

    TEST(FullFieldGrating_test_898){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 2, 5);
    }

    TEST(FullFieldGrating_test_899){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 2, 5);
    }

    TEST(FullFieldGrating_test_900){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 2, 6);
    }

    TEST(FullFieldGrating_test_901){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 2, 6);
    }

    TEST(FullFieldGrating_test_902){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 2, 7);
    }

    TEST(FullFieldGrating_test_903){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 2, 7);
    }

    TEST(FullFieldGrating_test_904){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 3, 1);
    }

    TEST(FullFieldGrating_test_905){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 3, 1);
    }

    TEST(FullFieldGrating_test_906){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 3, 2);
    }

    TEST(FullFieldGrating_test_907){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 3, 2);
    }

    TEST(FullFieldGrating_test_908){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 3, 3);
    }

    TEST(FullFieldGrating_test_909){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 3, 3);
    }

    TEST(FullFieldGrating_test_910){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 3, 5);
    }

    TEST(FullFieldGrating_test_911){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 3, 5);
    }

    TEST(FullFieldGrating_test_912){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 3, 6);
    }

    TEST(FullFieldGrating_test_913){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 3, 6);
    }

    TEST(FullFieldGrating_test_914){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 3, 7);
    }

    TEST(FullFieldGrating_test_915){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 3, 7);
    }

    TEST(FullFieldGrating_test_916){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 5, 1);
    }

    TEST(FullFieldGrating_test_917){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 5, 1);
    }

    TEST(FullFieldGrating_test_918){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 5, 2);
    }

    TEST(FullFieldGrating_test_919){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 5, 2);
    }

    TEST(FullFieldGrating_test_920){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 5, 3);
    }

    TEST(FullFieldGrating_test_921){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 5, 3);
    }

    TEST(FullFieldGrating_test_922){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 5, 5);
    }

    TEST(FullFieldGrating_test_923){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 5, 5);
    }

    TEST(FullFieldGrating_test_924){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 5, 6);
    }

    TEST(FullFieldGrating_test_925){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 5, 6);
    }

    TEST(FullFieldGrating_test_926){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 5, 7);
    }

    TEST(FullFieldGrating_test_927){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 5, 7);
    }

    TEST(FullFieldGrating_test_928){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 6, 1);
    }

    TEST(FullFieldGrating_test_929){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 6, 1);
    }

    TEST(FullFieldGrating_test_930){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 6, 2);
    }

    TEST(FullFieldGrating_test_931){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 6, 2);
    }

    TEST(FullFieldGrating_test_932){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 6, 3);
    }

    TEST(FullFieldGrating_test_933){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 6, 3);
    }

    TEST(FullFieldGrating_test_934){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 6, 5);
    }

    TEST(FullFieldGrating_test_935){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 6, 5);
    }

    TEST(FullFieldGrating_test_936){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 6, 6);
    }

    TEST(FullFieldGrating_test_937){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 6, 6);
    }

    TEST(FullFieldGrating_test_938){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 6, 7);
    }

    TEST(FullFieldGrating_test_939){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 6, 7);
    }

    TEST(FullFieldGrating_test_940){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 7, 1);
    }

    TEST(FullFieldGrating_test_941){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 7, 1);
    }

    TEST(FullFieldGrating_test_942){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 7, 2);
    }

    TEST(FullFieldGrating_test_943){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 7, 2);
    }

    TEST(FullFieldGrating_test_944){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 7, 3);
    }

    TEST(FullFieldGrating_test_945){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 7, 3);
    }

    TEST(FullFieldGrating_test_946){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 7, 5);
    }

    TEST(FullFieldGrating_test_947){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 7, 5);
    }

    TEST(FullFieldGrating_test_948){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 7, 6);
    }

    TEST(FullFieldGrating_test_949){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 7, 6);
    }

    TEST(FullFieldGrating_test_950){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 5, 7, 7);
    }

    TEST(FullFieldGrating_test_951){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 5, 7, 7);
    }

    TEST(FullFieldGrating_test_952){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 1, 1);
    }

    TEST(FullFieldGrating_test_953){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 1, 1);
    }

    TEST(FullFieldGrating_test_954){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 1, 2);
    }

    TEST(FullFieldGrating_test_955){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 1, 2);
    }

    TEST(FullFieldGrating_test_956){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 1, 3);
    }

    TEST(FullFieldGrating_test_957){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 1, 3);
    }

    TEST(FullFieldGrating_test_958){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 1, 5);
    }

    TEST(FullFieldGrating_test_959){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 1, 5);
    }

    TEST(FullFieldGrating_test_960){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 1, 6);
    }

    TEST(FullFieldGrating_test_961){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 1, 6);
    }

    TEST(FullFieldGrating_test_962){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 1, 7);
    }

    TEST(FullFieldGrating_test_963){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 1, 7);
    }

    TEST(FullFieldGrating_test_964){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 2, 1);
    }

    TEST(FullFieldGrating_test_965){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 2, 1);
    }

    TEST(FullFieldGrating_test_966){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 2, 2);
    }

    TEST(FullFieldGrating_test_967){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 2, 2);
    }

    TEST(FullFieldGrating_test_968){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 2, 3);
    }

    TEST(FullFieldGrating_test_969){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 2, 3);
    }

    TEST(FullFieldGrating_test_970){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 2, 5);
    }

    TEST(FullFieldGrating_test_971){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 2, 5);
    }

    TEST(FullFieldGrating_test_972){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 2, 6);
    }

    TEST(FullFieldGrating_test_973){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 2, 6);
    }

    TEST(FullFieldGrating_test_974){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 2, 7);
    }

    TEST(FullFieldGrating_test_975){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 2, 7);
    }

    TEST(FullFieldGrating_test_976){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 3, 1);
    }

    TEST(FullFieldGrating_test_977){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 3, 1);
    }

    TEST(FullFieldGrating_test_978){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 3, 2);
    }

    TEST(FullFieldGrating_test_979){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 3, 2);
    }

    TEST(FullFieldGrating_test_980){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 3, 3);
    }

    TEST(FullFieldGrating_test_981){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 3, 3);
    }

    TEST(FullFieldGrating_test_982){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 3, 5);
    }

    TEST(FullFieldGrating_test_983){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 3, 5);
    }

    TEST(FullFieldGrating_test_984){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 3, 6);
    }

    TEST(FullFieldGrating_test_985){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 3, 6);
    }

    TEST(FullFieldGrating_test_986){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 3, 7);
    }

    TEST(FullFieldGrating_test_987){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 3, 7);
    }

    TEST(FullFieldGrating_test_988){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 5, 1);
    }

    TEST(FullFieldGrating_test_989){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 5, 1);
    }

    TEST(FullFieldGrating_test_990){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 5, 2);
    }

    TEST(FullFieldGrating_test_991){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 5, 2);
    }

    TEST(FullFieldGrating_test_992){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 5, 3);
    }

    TEST(FullFieldGrating_test_993){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 5, 3);
    }

    TEST(FullFieldGrating_test_994){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 5, 5);
    }

    TEST(FullFieldGrating_test_995){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 5, 5);
    }

    TEST(FullFieldGrating_test_996){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 5, 6);
    }

    TEST(FullFieldGrating_test_997){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 5, 6);
    }

    TEST(FullFieldGrating_test_998){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 5, 7);
    }

    TEST(FullFieldGrating_test_999){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 5, 7);
    }

    TEST(FullFieldGrating_test_1000){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 6, 1);
    }

    TEST(FullFieldGrating_test_1001){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 6, 1);
    }

    TEST(FullFieldGrating_test_1002){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 6, 2);
    }

    TEST(FullFieldGrating_test_1003){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 6, 2);
    }

    TEST(FullFieldGrating_test_1004){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 6, 3);
    }

    TEST(FullFieldGrating_test_1005){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 6, 3);
    }

    TEST(FullFieldGrating_test_1006){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 6, 5);
    }

    TEST(FullFieldGrating_test_1007){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 6, 5);
    }

    TEST(FullFieldGrating_test_1008){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 6, 6);
    }

    TEST(FullFieldGrating_test_1009){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 6, 6);
    }

    TEST(FullFieldGrating_test_1010){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 6, 7);
    }

    TEST(FullFieldGrating_test_1011){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 6, 7);
    }

    TEST(FullFieldGrating_test_1012){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 7, 1);
    }

    TEST(FullFieldGrating_test_1013){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 7, 1);
    }

    TEST(FullFieldGrating_test_1014){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 7, 2);
    }

    TEST(FullFieldGrating_test_1015){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 7, 2);
    }

    TEST(FullFieldGrating_test_1016){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 7, 3);
    }

    TEST(FullFieldGrating_test_1017){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 7, 3);
    }

    TEST(FullFieldGrating_test_1018){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 7, 5);
    }

    TEST(FullFieldGrating_test_1019){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 7, 5);
    }

    TEST(FullFieldGrating_test_1020){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 7, 6);
    }

    TEST(FullFieldGrating_test_1021){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 7, 6);
    }

    TEST(FullFieldGrating_test_1022){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 6, 7, 7);
    }

    TEST(FullFieldGrating_test_1023){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 6, 7, 7);
    }

    TEST(FullFieldGrating_test_1024){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 1, 1);
    }

    TEST(FullFieldGrating_test_1025){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 1, 1);
    }

    TEST(FullFieldGrating_test_1026){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 1, 2);
    }

    TEST(FullFieldGrating_test_1027){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 1, 2);
    }

    TEST(FullFieldGrating_test_1028){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 1, 3);
    }

    TEST(FullFieldGrating_test_1029){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 1, 3);
    }

    TEST(FullFieldGrating_test_1030){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 1, 5);
    }

    TEST(FullFieldGrating_test_1031){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 1, 5);
    }

    TEST(FullFieldGrating_test_1032){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 1, 6);
    }

    TEST(FullFieldGrating_test_1033){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 1, 6);
    }

    TEST(FullFieldGrating_test_1034){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 1, 7);
    }

    TEST(FullFieldGrating_test_1035){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 1, 7);
    }

    TEST(FullFieldGrating_test_1036){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 2, 1);
    }

    TEST(FullFieldGrating_test_1037){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 2, 1);
    }

    TEST(FullFieldGrating_test_1038){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 2, 2);
    }

    TEST(FullFieldGrating_test_1039){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 2, 2);
    }

    TEST(FullFieldGrating_test_1040){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 2, 3);
    }

    TEST(FullFieldGrating_test_1041){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 2, 3);
    }

    TEST(FullFieldGrating_test_1042){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 2, 5);
    }

    TEST(FullFieldGrating_test_1043){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 2, 5);
    }

    TEST(FullFieldGrating_test_1044){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 2, 6);
    }

    TEST(FullFieldGrating_test_1045){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 2, 6);
    }

    TEST(FullFieldGrating_test_1046){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 2, 7);
    }

    TEST(FullFieldGrating_test_1047){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 2, 7);
    }

    TEST(FullFieldGrating_test_1048){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 3, 1);
    }

    TEST(FullFieldGrating_test_1049){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 3, 1);
    }

    TEST(FullFieldGrating_test_1050){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 3, 2);
    }

    TEST(FullFieldGrating_test_1051){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 3, 2);
    }

    TEST(FullFieldGrating_test_1052){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 3, 3);
    }

    TEST(FullFieldGrating_test_1053){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 3, 3);
    }

    TEST(FullFieldGrating_test_1054){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 3, 5);
    }

    TEST(FullFieldGrating_test_1055){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 3, 5);
    }

    TEST(FullFieldGrating_test_1056){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 3, 6);
    }

    TEST(FullFieldGrating_test_1057){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 3, 6);
    }

    TEST(FullFieldGrating_test_1058){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 3, 7);
    }

    TEST(FullFieldGrating_test_1059){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 3, 7);
    }

    TEST(FullFieldGrating_test_1060){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 5, 1);
    }

    TEST(FullFieldGrating_test_1061){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 5, 1);
    }

    TEST(FullFieldGrating_test_1062){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 5, 2);
    }

    TEST(FullFieldGrating_test_1063){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 5, 2);
    }

    TEST(FullFieldGrating_test_1064){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 5, 3);
    }

    TEST(FullFieldGrating_test_1065){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 5, 3);
    }

    TEST(FullFieldGrating_test_1066){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 5, 5);
    }

    TEST(FullFieldGrating_test_1067){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 5, 5);
    }

    TEST(FullFieldGrating_test_1068){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 5, 6);
    }

    TEST(FullFieldGrating_test_1069){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 5, 6);
    }

    TEST(FullFieldGrating_test_1070){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 5, 7);
    }

    TEST(FullFieldGrating_test_1071){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 5, 7);
    }

    TEST(FullFieldGrating_test_1072){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 6, 1);
    }

    TEST(FullFieldGrating_test_1073){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 6, 1);
    }

    TEST(FullFieldGrating_test_1074){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 6, 2);
    }

    TEST(FullFieldGrating_test_1075){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 6, 2);
    }

    TEST(FullFieldGrating_test_1076){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 6, 3);
    }

    TEST(FullFieldGrating_test_1077){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 6, 3);
    }

    TEST(FullFieldGrating_test_1078){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 6, 5);
    }

    TEST(FullFieldGrating_test_1079){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 6, 5);
    }

    TEST(FullFieldGrating_test_1080){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 6, 6);
    }

    TEST(FullFieldGrating_test_1081){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 6, 6);
    }

    TEST(FullFieldGrating_test_1082){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 6, 7);
    }

    TEST(FullFieldGrating_test_1083){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 6, 7);
    }

    TEST(FullFieldGrating_test_1084){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 7, 1);
    }

    TEST(FullFieldGrating_test_1085){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 7, 1);
    }

    TEST(FullFieldGrating_test_1086){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 7, 2);
    }

    TEST(FullFieldGrating_test_1087){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 7, 2);
    }

    TEST(FullFieldGrating_test_1088){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 7, 3);
    }

    TEST(FullFieldGrating_test_1089){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 7, 3);
    }

    TEST(FullFieldGrating_test_1090){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 7, 5);
    }

    TEST(FullFieldGrating_test_1091){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 7, 5);
    }

    TEST(FullFieldGrating_test_1092){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 7, 6);
    }

    TEST(FullFieldGrating_test_1093){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 7, 6);
    }

    TEST(FullFieldGrating_test_1094){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, -0.1, 7, 7, 7);
    }

    TEST(FullFieldGrating_test_1095){
         runStimulusFullFieldGratingTest(3, 3, 0.01, 0.01, 2.2, 7, 7, 7);
    }

    TEST(FullFieldGrating_test_1096){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 1, 1);
    }

    TEST(FullFieldGrating_test_1097){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 1, 1);
    }

    TEST(FullFieldGrating_test_1098){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 1, 2);
    }

    TEST(FullFieldGrating_test_1099){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 1, 2);
    }

    TEST(FullFieldGrating_test_1100){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 1, 3);
    }

    TEST(FullFieldGrating_test_1101){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 1, 3);
    }

    TEST(FullFieldGrating_test_1102){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 1, 5);
    }

    TEST(FullFieldGrating_test_1103){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 1, 5);
    }

    TEST(FullFieldGrating_test_1104){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 1, 6);
    }

    TEST(FullFieldGrating_test_1105){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 1, 6);
    }

    TEST(FullFieldGrating_test_1106){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 1, 7);
    }

    TEST(FullFieldGrating_test_1107){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 1, 7);
    }

    TEST(FullFieldGrating_test_1108){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 2, 1);
    }

    TEST(FullFieldGrating_test_1109){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 2, 1);
    }

    TEST(FullFieldGrating_test_1110){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 2, 2);
    }

    TEST(FullFieldGrating_test_1111){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 2, 2);
    }

    TEST(FullFieldGrating_test_1112){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 2, 3);
    }

    TEST(FullFieldGrating_test_1113){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 2, 3);
    }

    TEST(FullFieldGrating_test_1114){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 2, 5);
    }

    TEST(FullFieldGrating_test_1115){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 2, 5);
    }

    TEST(FullFieldGrating_test_1116){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 2, 6);
    }

    TEST(FullFieldGrating_test_1117){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 2, 6);
    }

    TEST(FullFieldGrating_test_1118){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 2, 7);
    }

    TEST(FullFieldGrating_test_1119){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 2, 7);
    }

    TEST(FullFieldGrating_test_1120){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 3, 1);
    }

    TEST(FullFieldGrating_test_1121){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 3, 1);
    }

    TEST(FullFieldGrating_test_1122){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 3, 2);
    }

    TEST(FullFieldGrating_test_1123){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 3, 2);
    }

    TEST(FullFieldGrating_test_1124){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 3, 3);
    }

    TEST(FullFieldGrating_test_1125){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 3, 3);
    }

    TEST(FullFieldGrating_test_1126){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 3, 5);
    }

    TEST(FullFieldGrating_test_1127){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 3, 5);
    }

    TEST(FullFieldGrating_test_1128){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 3, 6);
    }

    TEST(FullFieldGrating_test_1129){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 3, 6);
    }

    TEST(FullFieldGrating_test_1130){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 3, 7);
    }

    TEST(FullFieldGrating_test_1131){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 3, 7);
    }

    TEST(FullFieldGrating_test_1132){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 5, 1);
    }

    TEST(FullFieldGrating_test_1133){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 5, 1);
    }

    TEST(FullFieldGrating_test_1134){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 5, 2);
    }

    TEST(FullFieldGrating_test_1135){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 5, 2);
    }

    TEST(FullFieldGrating_test_1136){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 5, 3);
    }

    TEST(FullFieldGrating_test_1137){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 5, 3);
    }

    TEST(FullFieldGrating_test_1138){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 5, 5);
    }

    TEST(FullFieldGrating_test_1139){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 5, 5);
    }

    TEST(FullFieldGrating_test_1140){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 5, 6);
    }

    TEST(FullFieldGrating_test_1141){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 5, 6);
    }

    TEST(FullFieldGrating_test_1142){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 5, 7);
    }

    TEST(FullFieldGrating_test_1143){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 5, 7);
    }

    TEST(FullFieldGrating_test_1144){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 6, 1);
    }

    TEST(FullFieldGrating_test_1145){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 6, 1);
    }

    TEST(FullFieldGrating_test_1146){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 6, 2);
    }

    TEST(FullFieldGrating_test_1147){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 6, 2);
    }

    TEST(FullFieldGrating_test_1148){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 6, 3);
    }

    TEST(FullFieldGrating_test_1149){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 6, 3);
    }

    TEST(FullFieldGrating_test_1150){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 6, 5);
    }

    TEST(FullFieldGrating_test_1151){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 6, 5);
    }

    TEST(FullFieldGrating_test_1152){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 6, 6);
    }

    TEST(FullFieldGrating_test_1153){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 6, 6);
    }

    TEST(FullFieldGrating_test_1154){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 6, 7);
    }

    TEST(FullFieldGrating_test_1155){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 6, 7);
    }

    TEST(FullFieldGrating_test_1156){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 7, 1);
    }

    TEST(FullFieldGrating_test_1157){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 7, 1);
    }

    TEST(FullFieldGrating_test_1158){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 7, 2);
    }

    TEST(FullFieldGrating_test_1159){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 7, 2);
    }

    TEST(FullFieldGrating_test_1160){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 7, 3);
    }

    TEST(FullFieldGrating_test_1161){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 7, 3);
    }

    TEST(FullFieldGrating_test_1162){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 7, 5);
    }

    TEST(FullFieldGrating_test_1163){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 7, 5);
    }

    TEST(FullFieldGrating_test_1164){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 7, 6);
    }

    TEST(FullFieldGrating_test_1165){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 7, 6);
    }

    TEST(FullFieldGrating_test_1166){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 0, 7, 7);
    }

    TEST(FullFieldGrating_test_1167){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 0, 7, 7);
    }

    TEST(FullFieldGrating_test_1168){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 1, 1);
    }

    TEST(FullFieldGrating_test_1169){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 1, 1);
    }

    TEST(FullFieldGrating_test_1170){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 1, 2);
    }

    TEST(FullFieldGrating_test_1171){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 1, 2);
    }

    TEST(FullFieldGrating_test_1172){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 1, 3);
    }

    TEST(FullFieldGrating_test_1173){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 1, 3);
    }

    TEST(FullFieldGrating_test_1174){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 1, 5);
    }

    TEST(FullFieldGrating_test_1175){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 1, 5);
    }

    TEST(FullFieldGrating_test_1176){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 1, 6);
    }

    TEST(FullFieldGrating_test_1177){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 1, 6);
    }

    TEST(FullFieldGrating_test_1178){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 1, 7);
    }

    TEST(FullFieldGrating_test_1179){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 1, 7);
    }

    TEST(FullFieldGrating_test_1180){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 2, 1);
    }

    TEST(FullFieldGrating_test_1181){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 2, 1);
    }

    TEST(FullFieldGrating_test_1182){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 2, 2);
    }

    TEST(FullFieldGrating_test_1183){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 2, 2);
    }

    TEST(FullFieldGrating_test_1184){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 2, 3);
    }

    TEST(FullFieldGrating_test_1185){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 2, 3);
    }

    TEST(FullFieldGrating_test_1186){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 2, 5);
    }

    TEST(FullFieldGrating_test_1187){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 2, 5);
    }

    TEST(FullFieldGrating_test_1188){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 2, 6);
    }

    TEST(FullFieldGrating_test_1189){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 2, 6);
    }

    TEST(FullFieldGrating_test_1190){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 2, 7);
    }

    TEST(FullFieldGrating_test_1191){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 2, 7);
    }

    TEST(FullFieldGrating_test_1192){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 3, 1);
    }

    TEST(FullFieldGrating_test_1193){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 3, 1);
    }

    TEST(FullFieldGrating_test_1194){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 3, 2);
    }

    TEST(FullFieldGrating_test_1195){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 3, 2);
    }

    TEST(FullFieldGrating_test_1196){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 3, 3);
    }

    TEST(FullFieldGrating_test_1197){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 3, 3);
    }

    TEST(FullFieldGrating_test_1198){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 3, 5);
    }

    TEST(FullFieldGrating_test_1199){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 3, 5);
    }

    TEST(FullFieldGrating_test_1200){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 3, 6);
    }

    TEST(FullFieldGrating_test_1201){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 3, 6);
    }

    TEST(FullFieldGrating_test_1202){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 3, 7);
    }

    TEST(FullFieldGrating_test_1203){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 3, 7);
    }

    TEST(FullFieldGrating_test_1204){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 5, 1);
    }

    TEST(FullFieldGrating_test_1205){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 5, 1);
    }

    TEST(FullFieldGrating_test_1206){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 5, 2);
    }

    TEST(FullFieldGrating_test_1207){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 5, 2);
    }

    TEST(FullFieldGrating_test_1208){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 5, 3);
    }

    TEST(FullFieldGrating_test_1209){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 5, 3);
    }

    TEST(FullFieldGrating_test_1210){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 5, 5);
    }

    TEST(FullFieldGrating_test_1211){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 5, 5);
    }

    TEST(FullFieldGrating_test_1212){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 5, 6);
    }

    TEST(FullFieldGrating_test_1213){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 5, 6);
    }

    TEST(FullFieldGrating_test_1214){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 5, 7);
    }

    TEST(FullFieldGrating_test_1215){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 5, 7);
    }

    TEST(FullFieldGrating_test_1216){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 6, 1);
    }

    TEST(FullFieldGrating_test_1217){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 6, 1);
    }

    TEST(FullFieldGrating_test_1218){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 6, 2);
    }

    TEST(FullFieldGrating_test_1219){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 6, 2);
    }

    TEST(FullFieldGrating_test_1220){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 6, 3);
    }

    TEST(FullFieldGrating_test_1221){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 6, 3);
    }

    TEST(FullFieldGrating_test_1222){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 6, 5);
    }

    TEST(FullFieldGrating_test_1223){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 6, 5);
    }

    TEST(FullFieldGrating_test_1224){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 6, 6);
    }

    TEST(FullFieldGrating_test_1225){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 6, 6);
    }

    TEST(FullFieldGrating_test_1226){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 6, 7);
    }

    TEST(FullFieldGrating_test_1227){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 6, 7);
    }

    TEST(FullFieldGrating_test_1228){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 7, 1);
    }

    TEST(FullFieldGrating_test_1229){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 7, 1);
    }

    TEST(FullFieldGrating_test_1230){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 7, 2);
    }

    TEST(FullFieldGrating_test_1231){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 7, 2);
    }

    TEST(FullFieldGrating_test_1232){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 7, 3);
    }

    TEST(FullFieldGrating_test_1233){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 7, 3);
    }

    TEST(FullFieldGrating_test_1234){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 7, 5);
    }

    TEST(FullFieldGrating_test_1235){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 7, 5);
    }

    TEST(FullFieldGrating_test_1236){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 7, 6);
    }

    TEST(FullFieldGrating_test_1237){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 7, 6);
    }

    TEST(FullFieldGrating_test_1238){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 1, 7, 7);
    }

    TEST(FullFieldGrating_test_1239){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 1, 7, 7);
    }

    TEST(FullFieldGrating_test_1240){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 1, 1);
    }

    TEST(FullFieldGrating_test_1241){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 1, 1);
    }

    TEST(FullFieldGrating_test_1242){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 1, 2);
    }

    TEST(FullFieldGrating_test_1243){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 1, 2);
    }

    TEST(FullFieldGrating_test_1244){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 1, 3);
    }

    TEST(FullFieldGrating_test_1245){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 1, 3);
    }

    TEST(FullFieldGrating_test_1246){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 1, 5);
    }

    TEST(FullFieldGrating_test_1247){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 1, 5);
    }

    TEST(FullFieldGrating_test_1248){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 1, 6);
    }

    TEST(FullFieldGrating_test_1249){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 1, 6);
    }

    TEST(FullFieldGrating_test_1250){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 1, 7);
    }

    TEST(FullFieldGrating_test_1251){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 1, 7);
    }

    TEST(FullFieldGrating_test_1252){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 2, 1);
    }

    TEST(FullFieldGrating_test_1253){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 2, 1);
    }

    TEST(FullFieldGrating_test_1254){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 2, 2);
    }

    TEST(FullFieldGrating_test_1255){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 2, 2);
    }

    TEST(FullFieldGrating_test_1256){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 2, 3);
    }

    TEST(FullFieldGrating_test_1257){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 2, 3);
    }

    TEST(FullFieldGrating_test_1258){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 2, 5);
    }

    TEST(FullFieldGrating_test_1259){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 2, 5);
    }

    TEST(FullFieldGrating_test_1260){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 2, 6);
    }

    TEST(FullFieldGrating_test_1261){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 2, 6);
    }

    TEST(FullFieldGrating_test_1262){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 2, 7);
    }

    TEST(FullFieldGrating_test_1263){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 2, 7);
    }

    TEST(FullFieldGrating_test_1264){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 3, 1);
    }

    TEST(FullFieldGrating_test_1265){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 3, 1);
    }

    TEST(FullFieldGrating_test_1266){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 3, 2);
    }

    TEST(FullFieldGrating_test_1267){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 3, 2);
    }

    TEST(FullFieldGrating_test_1268){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 3, 3);
    }

    TEST(FullFieldGrating_test_1269){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 3, 3);
    }

    TEST(FullFieldGrating_test_1270){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 3, 5);
    }

    TEST(FullFieldGrating_test_1271){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 3, 5);
    }

    TEST(FullFieldGrating_test_1272){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 3, 6);
    }

    TEST(FullFieldGrating_test_1273){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 3, 6);
    }

    TEST(FullFieldGrating_test_1274){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 3, 7);
    }

    TEST(FullFieldGrating_test_1275){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 3, 7);
    }

    TEST(FullFieldGrating_test_1276){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 5, 1);
    }

    TEST(FullFieldGrating_test_1277){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 5, 1);
    }

    TEST(FullFieldGrating_test_1278){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 5, 2);
    }

    TEST(FullFieldGrating_test_1279){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 5, 2);
    }

    TEST(FullFieldGrating_test_1280){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 5, 3);
    }

    TEST(FullFieldGrating_test_1281){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 5, 3);
    }

    TEST(FullFieldGrating_test_1282){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 5, 5);
    }

    TEST(FullFieldGrating_test_1283){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 5, 5);
    }

    TEST(FullFieldGrating_test_1284){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 5, 6);
    }

    TEST(FullFieldGrating_test_1285){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 5, 6);
    }

    TEST(FullFieldGrating_test_1286){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 5, 7);
    }

    TEST(FullFieldGrating_test_1287){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 5, 7);
    }

    TEST(FullFieldGrating_test_1288){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 6, 1);
    }

    TEST(FullFieldGrating_test_1289){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 6, 1);
    }

    TEST(FullFieldGrating_test_1290){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 6, 2);
    }

    TEST(FullFieldGrating_test_1291){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 6, 2);
    }

    TEST(FullFieldGrating_test_1292){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 6, 3);
    }

    TEST(FullFieldGrating_test_1293){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 6, 3);
    }

    TEST(FullFieldGrating_test_1294){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 6, 5);
    }

    TEST(FullFieldGrating_test_1295){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 6, 5);
    }

    TEST(FullFieldGrating_test_1296){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 6, 6);
    }

    TEST(FullFieldGrating_test_1297){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 6, 6);
    }

    TEST(FullFieldGrating_test_1298){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 6, 7);
    }

    TEST(FullFieldGrating_test_1299){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 6, 7);
    }

    TEST(FullFieldGrating_test_1300){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 7, 1);
    }

    TEST(FullFieldGrating_test_1301){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 7, 1);
    }

    TEST(FullFieldGrating_test_1302){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 7, 2);
    }

    TEST(FullFieldGrating_test_1303){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 7, 2);
    }

    TEST(FullFieldGrating_test_1304){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 7, 3);
    }

    TEST(FullFieldGrating_test_1305){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 7, 3);
    }

    TEST(FullFieldGrating_test_1306){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 7, 5);
    }

    TEST(FullFieldGrating_test_1307){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 7, 5);
    }

    TEST(FullFieldGrating_test_1308){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 7, 6);
    }

    TEST(FullFieldGrating_test_1309){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 7, 6);
    }

    TEST(FullFieldGrating_test_1310){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 2, 7, 7);
    }

    TEST(FullFieldGrating_test_1311){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 2, 7, 7);
    }

    TEST(FullFieldGrating_test_1312){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 1, 1);
    }

    TEST(FullFieldGrating_test_1313){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 1, 1);
    }

    TEST(FullFieldGrating_test_1314){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 1, 2);
    }

    TEST(FullFieldGrating_test_1315){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 1, 2);
    }

    TEST(FullFieldGrating_test_1316){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 1, 3);
    }

    TEST(FullFieldGrating_test_1317){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 1, 3);
    }

    TEST(FullFieldGrating_test_1318){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 1, 5);
    }

    TEST(FullFieldGrating_test_1319){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 1, 5);
    }

    TEST(FullFieldGrating_test_1320){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 1, 6);
    }

    TEST(FullFieldGrating_test_1321){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 1, 6);
    }

    TEST(FullFieldGrating_test_1322){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 1, 7);
    }

    TEST(FullFieldGrating_test_1323){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 1, 7);
    }

    TEST(FullFieldGrating_test_1324){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 2, 1);
    }

    TEST(FullFieldGrating_test_1325){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 2, 1);
    }

    TEST(FullFieldGrating_test_1326){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 2, 2);
    }

    TEST(FullFieldGrating_test_1327){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 2, 2);
    }

    TEST(FullFieldGrating_test_1328){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 2, 3);
    }

    TEST(FullFieldGrating_test_1329){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 2, 3);
    }

    TEST(FullFieldGrating_test_1330){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 2, 5);
    }

    TEST(FullFieldGrating_test_1331){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 2, 5);
    }

    TEST(FullFieldGrating_test_1332){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 2, 6);
    }

    TEST(FullFieldGrating_test_1333){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 2, 6);
    }

    TEST(FullFieldGrating_test_1334){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 2, 7);
    }

    TEST(FullFieldGrating_test_1335){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 2, 7);
    }

    TEST(FullFieldGrating_test_1336){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 3, 1);
    }

    TEST(FullFieldGrating_test_1337){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 3, 1);
    }

    TEST(FullFieldGrating_test_1338){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 3, 2);
    }

    TEST(FullFieldGrating_test_1339){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 3, 2);
    }

    TEST(FullFieldGrating_test_1340){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 3, 3);
    }

    TEST(FullFieldGrating_test_1341){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 3, 3);
    }

    TEST(FullFieldGrating_test_1342){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 3, 5);
    }

    TEST(FullFieldGrating_test_1343){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 3, 5);
    }

    TEST(FullFieldGrating_test_1344){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 3, 6);
    }

    TEST(FullFieldGrating_test_1345){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 3, 6);
    }

    TEST(FullFieldGrating_test_1346){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 3, 7);
    }

    TEST(FullFieldGrating_test_1347){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 3, 7);
    }

    TEST(FullFieldGrating_test_1348){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 5, 1);
    }

    TEST(FullFieldGrating_test_1349){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 5, 1);
    }

    TEST(FullFieldGrating_test_1350){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 5, 2);
    }

    TEST(FullFieldGrating_test_1351){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 5, 2);
    }

    TEST(FullFieldGrating_test_1352){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 5, 3);
    }

    TEST(FullFieldGrating_test_1353){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 5, 3);
    }

    TEST(FullFieldGrating_test_1354){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 5, 5);
    }

    TEST(FullFieldGrating_test_1355){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 5, 5);
    }

    TEST(FullFieldGrating_test_1356){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 5, 6);
    }

    TEST(FullFieldGrating_test_1357){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 5, 6);
    }

    TEST(FullFieldGrating_test_1358){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 5, 7);
    }

    TEST(FullFieldGrating_test_1359){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 5, 7);
    }

    TEST(FullFieldGrating_test_1360){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 6, 1);
    }

    TEST(FullFieldGrating_test_1361){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 6, 1);
    }

    TEST(FullFieldGrating_test_1362){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 6, 2);
    }

    TEST(FullFieldGrating_test_1363){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 6, 2);
    }

    TEST(FullFieldGrating_test_1364){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 6, 3);
    }

    TEST(FullFieldGrating_test_1365){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 6, 3);
    }

    TEST(FullFieldGrating_test_1366){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 6, 5);
    }

    TEST(FullFieldGrating_test_1367){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 6, 5);
    }

    TEST(FullFieldGrating_test_1368){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 6, 6);
    }

    TEST(FullFieldGrating_test_1369){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 6, 6);
    }

    TEST(FullFieldGrating_test_1370){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 6, 7);
    }

    TEST(FullFieldGrating_test_1371){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 6, 7);
    }

    TEST(FullFieldGrating_test_1372){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 7, 1);
    }

    TEST(FullFieldGrating_test_1373){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 7, 1);
    }

    TEST(FullFieldGrating_test_1374){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 7, 2);
    }

    TEST(FullFieldGrating_test_1375){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 7, 2);
    }

    TEST(FullFieldGrating_test_1376){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 7, 3);
    }

    TEST(FullFieldGrating_test_1377){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 7, 3);
    }

    TEST(FullFieldGrating_test_1378){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 7, 5);
    }

    TEST(FullFieldGrating_test_1379){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 7, 5);
    }

    TEST(FullFieldGrating_test_1380){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 7, 6);
    }

    TEST(FullFieldGrating_test_1381){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 7, 6);
    }

    TEST(FullFieldGrating_test_1382){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 3, 7, 7);
    }

    TEST(FullFieldGrating_test_1383){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 3, 7, 7);
    }

    TEST(FullFieldGrating_test_1384){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 1, 1);
    }

    TEST(FullFieldGrating_test_1385){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 1, 1);
    }

    TEST(FullFieldGrating_test_1386){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 1, 2);
    }

    TEST(FullFieldGrating_test_1387){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 1, 2);
    }

    TEST(FullFieldGrating_test_1388){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 1, 3);
    }

    TEST(FullFieldGrating_test_1389){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 1, 3);
    }

    TEST(FullFieldGrating_test_1390){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 1, 5);
    }

    TEST(FullFieldGrating_test_1391){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 1, 5);
    }

    TEST(FullFieldGrating_test_1392){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 1, 6);
    }

    TEST(FullFieldGrating_test_1393){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 1, 6);
    }

    TEST(FullFieldGrating_test_1394){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 1, 7);
    }

    TEST(FullFieldGrating_test_1395){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 1, 7);
    }

    TEST(FullFieldGrating_test_1396){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 2, 1);
    }

    TEST(FullFieldGrating_test_1397){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 2, 1);
    }

    TEST(FullFieldGrating_test_1398){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 2, 2);
    }

    TEST(FullFieldGrating_test_1399){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 2, 2);
    }

    TEST(FullFieldGrating_test_1400){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 2, 3);
    }

    TEST(FullFieldGrating_test_1401){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 2, 3);
    }

    TEST(FullFieldGrating_test_1402){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 2, 5);
    }

    TEST(FullFieldGrating_test_1403){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 2, 5);
    }

    TEST(FullFieldGrating_test_1404){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 2, 6);
    }

    TEST(FullFieldGrating_test_1405){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 2, 6);
    }

    TEST(FullFieldGrating_test_1406){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 2, 7);
    }

    TEST(FullFieldGrating_test_1407){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 2, 7);
    }

    TEST(FullFieldGrating_test_1408){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 3, 1);
    }

    TEST(FullFieldGrating_test_1409){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 3, 1);
    }

    TEST(FullFieldGrating_test_1410){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 3, 2);
    }

    TEST(FullFieldGrating_test_1411){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 3, 2);
    }

    TEST(FullFieldGrating_test_1412){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 3, 3);
    }

    TEST(FullFieldGrating_test_1413){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 3, 3);
    }

    TEST(FullFieldGrating_test_1414){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 3, 5);
    }

    TEST(FullFieldGrating_test_1415){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 3, 5);
    }

    TEST(FullFieldGrating_test_1416){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 3, 6);
    }

    TEST(FullFieldGrating_test_1417){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 3, 6);
    }

    TEST(FullFieldGrating_test_1418){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 3, 7);
    }

    TEST(FullFieldGrating_test_1419){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 3, 7);
    }

    TEST(FullFieldGrating_test_1420){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 5, 1);
    }

    TEST(FullFieldGrating_test_1421){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 5, 1);
    }

    TEST(FullFieldGrating_test_1422){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 5, 2);
    }

    TEST(FullFieldGrating_test_1423){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 5, 2);
    }

    TEST(FullFieldGrating_test_1424){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 5, 3);
    }

    TEST(FullFieldGrating_test_1425){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 5, 3);
    }

    TEST(FullFieldGrating_test_1426){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 5, 5);
    }

    TEST(FullFieldGrating_test_1427){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 5, 5);
    }

    TEST(FullFieldGrating_test_1428){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 5, 6);
    }

    TEST(FullFieldGrating_test_1429){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 5, 6);
    }

    TEST(FullFieldGrating_test_1430){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 5, 7);
    }

    TEST(FullFieldGrating_test_1431){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 5, 7);
    }

    TEST(FullFieldGrating_test_1432){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 6, 1);
    }

    TEST(FullFieldGrating_test_1433){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 6, 1);
    }

    TEST(FullFieldGrating_test_1434){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 6, 2);
    }

    TEST(FullFieldGrating_test_1435){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 6, 2);
    }

    TEST(FullFieldGrating_test_1436){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 6, 3);
    }

    TEST(FullFieldGrating_test_1437){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 6, 3);
    }

    TEST(FullFieldGrating_test_1438){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 6, 5);
    }

    TEST(FullFieldGrating_test_1439){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 6, 5);
    }

    TEST(FullFieldGrating_test_1440){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 6, 6);
    }

    TEST(FullFieldGrating_test_1441){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 6, 6);
    }

    TEST(FullFieldGrating_test_1442){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 6, 7);
    }

    TEST(FullFieldGrating_test_1443){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 6, 7);
    }

    TEST(FullFieldGrating_test_1444){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 7, 1);
    }

    TEST(FullFieldGrating_test_1445){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 7, 1);
    }

    TEST(FullFieldGrating_test_1446){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 7, 2);
    }

    TEST(FullFieldGrating_test_1447){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 7, 2);
    }

    TEST(FullFieldGrating_test_1448){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 7, 3);
    }

    TEST(FullFieldGrating_test_1449){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 7, 3);
    }

    TEST(FullFieldGrating_test_1450){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 7, 5);
    }

    TEST(FullFieldGrating_test_1451){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 7, 5);
    }

    TEST(FullFieldGrating_test_1452){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 7, 6);
    }

    TEST(FullFieldGrating_test_1453){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 7, 6);
    }

    TEST(FullFieldGrating_test_1454){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 5, 7, 7);
    }

    TEST(FullFieldGrating_test_1455){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 5, 7, 7);
    }

    TEST(FullFieldGrating_test_1456){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 1, 1);
    }

    TEST(FullFieldGrating_test_1457){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 1, 1);
    }

    TEST(FullFieldGrating_test_1458){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 1, 2);
    }

    TEST(FullFieldGrating_test_1459){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 1, 2);
    }

    TEST(FullFieldGrating_test_1460){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 1, 3);
    }

    TEST(FullFieldGrating_test_1461){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 1, 3);
    }

    TEST(FullFieldGrating_test_1462){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 1, 5);
    }

    TEST(FullFieldGrating_test_1463){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 1, 5);
    }

    TEST(FullFieldGrating_test_1464){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 1, 6);
    }

    TEST(FullFieldGrating_test_1465){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 1, 6);
    }

    TEST(FullFieldGrating_test_1466){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 1, 7);
    }

    TEST(FullFieldGrating_test_1467){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 1, 7);
    }

    TEST(FullFieldGrating_test_1468){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 2, 1);
    }

    TEST(FullFieldGrating_test_1469){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 2, 1);
    }

    TEST(FullFieldGrating_test_1470){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 2, 2);
    }

    TEST(FullFieldGrating_test_1471){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 2, 2);
    }

    TEST(FullFieldGrating_test_1472){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 2, 3);
    }

    TEST(FullFieldGrating_test_1473){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 2, 3);
    }

    TEST(FullFieldGrating_test_1474){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 2, 5);
    }

    TEST(FullFieldGrating_test_1475){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 2, 5);
    }

    TEST(FullFieldGrating_test_1476){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 2, 6);
    }

    TEST(FullFieldGrating_test_1477){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 2, 6);
    }

    TEST(FullFieldGrating_test_1478){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 2, 7);
    }

    TEST(FullFieldGrating_test_1479){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 2, 7);
    }

    TEST(FullFieldGrating_test_1480){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 3, 1);
    }

    TEST(FullFieldGrating_test_1481){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 3, 1);
    }

    TEST(FullFieldGrating_test_1482){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 3, 2);
    }

    TEST(FullFieldGrating_test_1483){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 3, 2);
    }

    TEST(FullFieldGrating_test_1484){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 3, 3);
    }

    TEST(FullFieldGrating_test_1485){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 3, 3);
    }

    TEST(FullFieldGrating_test_1486){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 3, 5);
    }

    TEST(FullFieldGrating_test_1487){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 3, 5);
    }

    TEST(FullFieldGrating_test_1488){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 3, 6);
    }

    TEST(FullFieldGrating_test_1489){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 3, 6);
    }

    TEST(FullFieldGrating_test_1490){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 3, 7);
    }

    TEST(FullFieldGrating_test_1491){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 3, 7);
    }

    TEST(FullFieldGrating_test_1492){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 5, 1);
    }

    TEST(FullFieldGrating_test_1493){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 5, 1);
    }

    TEST(FullFieldGrating_test_1494){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 5, 2);
    }

    TEST(FullFieldGrating_test_1495){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 5, 2);
    }

    TEST(FullFieldGrating_test_1496){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 5, 3);
    }

    TEST(FullFieldGrating_test_1497){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 5, 3);
    }

    TEST(FullFieldGrating_test_1498){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 5, 5);
    }

    TEST(FullFieldGrating_test_1499){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 5, 5);
    }

    TEST(FullFieldGrating_test_1500){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 5, 6);
    }

    TEST(FullFieldGrating_test_1501){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 5, 6);
    }

    TEST(FullFieldGrating_test_1502){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 5, 7);
    }

    TEST(FullFieldGrating_test_1503){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 5, 7);
    }

    TEST(FullFieldGrating_test_1504){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 6, 1);
    }

    TEST(FullFieldGrating_test_1505){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 6, 1);
    }

    TEST(FullFieldGrating_test_1506){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 6, 2);
    }

    TEST(FullFieldGrating_test_1507){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 6, 2);
    }

    TEST(FullFieldGrating_test_1508){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 6, 3);
    }

    TEST(FullFieldGrating_test_1509){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 6, 3);
    }

    TEST(FullFieldGrating_test_1510){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 6, 5);
    }

    TEST(FullFieldGrating_test_1511){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 6, 5);
    }

    TEST(FullFieldGrating_test_1512){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 6, 6);
    }

    TEST(FullFieldGrating_test_1513){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 6, 6);
    }

    TEST(FullFieldGrating_test_1514){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 6, 7);
    }

    TEST(FullFieldGrating_test_1515){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 6, 7);
    }

    TEST(FullFieldGrating_test_1516){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 7, 1);
    }

    TEST(FullFieldGrating_test_1517){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 7, 1);
    }

    TEST(FullFieldGrating_test_1518){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 7, 2);
    }

    TEST(FullFieldGrating_test_1519){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 7, 2);
    }

    TEST(FullFieldGrating_test_1520){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 7, 3);
    }

    TEST(FullFieldGrating_test_1521){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 7, 3);
    }

    TEST(FullFieldGrating_test_1522){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 7, 5);
    }

    TEST(FullFieldGrating_test_1523){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 7, 5);
    }

    TEST(FullFieldGrating_test_1524){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 7, 6);
    }

    TEST(FullFieldGrating_test_1525){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 7, 6);
    }

    TEST(FullFieldGrating_test_1526){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 6, 7, 7);
    }

    TEST(FullFieldGrating_test_1527){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 6, 7, 7);
    }

    TEST(FullFieldGrating_test_1528){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 1, 1);
    }

    TEST(FullFieldGrating_test_1529){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 1, 1);
    }

    TEST(FullFieldGrating_test_1530){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 1, 2);
    }

    TEST(FullFieldGrating_test_1531){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 1, 2);
    }

    TEST(FullFieldGrating_test_1532){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 1, 3);
    }

    TEST(FullFieldGrating_test_1533){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 1, 3);
    }

    TEST(FullFieldGrating_test_1534){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 1, 5);
    }

    TEST(FullFieldGrating_test_1535){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 1, 5);
    }

    TEST(FullFieldGrating_test_1536){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 1, 6);
    }

    TEST(FullFieldGrating_test_1537){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 1, 6);
    }

    TEST(FullFieldGrating_test_1538){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 1, 7);
    }

    TEST(FullFieldGrating_test_1539){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 1, 7);
    }

    TEST(FullFieldGrating_test_1540){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 2, 1);
    }

    TEST(FullFieldGrating_test_1541){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 2, 1);
    }

    TEST(FullFieldGrating_test_1542){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 2, 2);
    }

    TEST(FullFieldGrating_test_1543){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 2, 2);
    }

    TEST(FullFieldGrating_test_1544){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 2, 3);
    }

    TEST(FullFieldGrating_test_1545){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 2, 3);
    }

    TEST(FullFieldGrating_test_1546){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 2, 5);
    }

    TEST(FullFieldGrating_test_1547){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 2, 5);
    }

    TEST(FullFieldGrating_test_1548){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 2, 6);
    }

    TEST(FullFieldGrating_test_1549){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 2, 6);
    }

    TEST(FullFieldGrating_test_1550){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 2, 7);
    }

    TEST(FullFieldGrating_test_1551){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 2, 7);
    }

    TEST(FullFieldGrating_test_1552){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 3, 1);
    }

    TEST(FullFieldGrating_test_1553){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 3, 1);
    }

    TEST(FullFieldGrating_test_1554){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 3, 2);
    }

    TEST(FullFieldGrating_test_1555){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 3, 2);
    }

    TEST(FullFieldGrating_test_1556){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 3, 3);
    }

    TEST(FullFieldGrating_test_1557){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 3, 3);
    }

    TEST(FullFieldGrating_test_1558){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 3, 5);
    }

    TEST(FullFieldGrating_test_1559){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 3, 5);
    }

    TEST(FullFieldGrating_test_1560){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 3, 6);
    }

    TEST(FullFieldGrating_test_1561){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 3, 6);
    }

    TEST(FullFieldGrating_test_1562){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 3, 7);
    }

    TEST(FullFieldGrating_test_1563){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 3, 7);
    }

    TEST(FullFieldGrating_test_1564){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 5, 1);
    }

    TEST(FullFieldGrating_test_1565){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 5, 1);
    }

    TEST(FullFieldGrating_test_1566){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 5, 2);
    }

    TEST(FullFieldGrating_test_1567){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 5, 2);
    }

    TEST(FullFieldGrating_test_1568){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 5, 3);
    }

    TEST(FullFieldGrating_test_1569){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 5, 3);
    }

    TEST(FullFieldGrating_test_1570){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 5, 5);
    }

    TEST(FullFieldGrating_test_1571){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 5, 5);
    }

    TEST(FullFieldGrating_test_1572){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 5, 6);
    }

    TEST(FullFieldGrating_test_1573){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 5, 6);
    }

    TEST(FullFieldGrating_test_1574){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 5, 7);
    }

    TEST(FullFieldGrating_test_1575){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 5, 7);
    }

    TEST(FullFieldGrating_test_1576){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 6, 1);
    }

    TEST(FullFieldGrating_test_1577){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 6, 1);
    }

    TEST(FullFieldGrating_test_1578){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 6, 2);
    }

    TEST(FullFieldGrating_test_1579){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 6, 2);
    }

    TEST(FullFieldGrating_test_1580){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 6, 3);
    }

    TEST(FullFieldGrating_test_1581){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 6, 3);
    }

    TEST(FullFieldGrating_test_1582){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 6, 5);
    }

    TEST(FullFieldGrating_test_1583){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 6, 5);
    }

    TEST(FullFieldGrating_test_1584){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 6, 6);
    }

    TEST(FullFieldGrating_test_1585){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 6, 6);
    }

    TEST(FullFieldGrating_test_1586){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 6, 7);
    }

    TEST(FullFieldGrating_test_1587){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 6, 7);
    }

    TEST(FullFieldGrating_test_1588){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 7, 1);
    }

    TEST(FullFieldGrating_test_1589){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 7, 1);
    }

    TEST(FullFieldGrating_test_1590){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 7, 2);
    }

    TEST(FullFieldGrating_test_1591){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 7, 2);
    }

    TEST(FullFieldGrating_test_1592){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 7, 3);
    }

    TEST(FullFieldGrating_test_1593){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 7, 3);
    }

    TEST(FullFieldGrating_test_1594){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 7, 5);
    }

    TEST(FullFieldGrating_test_1595){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 7, 5);
    }

    TEST(FullFieldGrating_test_1596){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 7, 6);
    }

    TEST(FullFieldGrating_test_1597){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 7, 6);
    }

    TEST(FullFieldGrating_test_1598){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, -0.1, 7, 7, 7);
    }

    TEST(FullFieldGrating_test_1599){
         runStimulusFullFieldGratingTest(3, 3, 0.5, 0.01, 2.2, 7, 7, 7);
    }
}


























