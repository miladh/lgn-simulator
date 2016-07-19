/**********************************************************************
 *  Test: ganglion cell response for patch grating, in the case where
 *        the stimulation spot and receptive-field are concentric.
 *
 *  Analytic source: Einevoll et. al (2005)
 *
 * ********************************************************************/


#include <unittest++/UnitTest++.h>
#include <lgnSimulator.h>

using namespace lgnSimulator;



double X(double y, double z, unsigned int n)
{
    double res = 0;
    double expTerm = exp(-z*z*0.25)/ (4*y*y);

    for(int i=0; i < int(n); i++){
        unsigned int n_fac = Special::factorial(i);
        double term = pow(z * 0.5, 2*i);
        double hypGeom = Special::confluentHypergeometric(i+1, 2, -1./(4*y*y));
        res += n_fac * term * hypGeom;
    }

    return res * expTerm;
}

double dogResponse(double C, double maskSize, double kd,
                   double a, double b, double c, int n)
{
    double term1 = X(a/maskSize, a*kd, n);
    double term2 = c * X(b/maskSize, b*kd, n);

    return C*(term1 - term2);

}


void runSystemTest_G_pg(int nt, double dt, int ns, double ds,
                          double C, double maskSize, int kId,
                          double a, double b, double c,
                          double weight, unsigned int n)
{

    //integrator
    Integrator integrator(nt, dt, ns, ds);

    //stimulus
    double kd = integrator.spatialFreqVec()[kId];
    CircleMaskGrating stim(&integrator, kd, 0, 0, C, maskSize);
    stim.computeFourierTransform();

    //ganglion cell
    DOG Ws(a, b, c);
    TemporalDelta Wt(0, dt);
    SeparableKernel W(weight, &Ws, &Wt);

    GanglionCell ganglion(&integrator, W);

    //Compute analytic:
    double Rg_ex = dogResponse(stim.contrast(), stim.maskSize(), stim.spatialFreq(),
                               a, b, c, n);

    //Compute numerical
    ganglion.computeResponse(&stim);
    double Rgc = ganglion.response()(integrator.nPointsSpatial()/2,
                                     integrator.nPointsSpatial()/2,
                                     0);


    //Test
    cout << kd << endl;
    cout << Rg_ex - Rgc << endl;
    CHECK_CLOSE(Rg_ex, Rgc, 1e-3);
}




SUITE(system){


    //Grating---------------------------------------------
    TEST(runTest_G_pg_8){
        double d = 1.5;
        double kId = 6;
        runSystemTest_G_pg(1, 1.0, 10, 0.1,
                           1.0, d, kId,
                           0.62, 1.26, 0.85,
                           1, 4);
    }

    TEST(runTest_G_pg_9){
        double d = 5.5;
        double kId = 6;
        runSystemTest_G_pg(1, 1.0, 10, 0.1,
                           1.0, d, kId,
                           0.62, 1.26, 0.85,
                           1, 4);
    }


    TEST(runTest_G_pg_10){
        double d = 5.5;
        double kId = 90;
        runSystemTest_G_pg(1, 1.0, 9, 0.1,
                           1.0, d, kId,
                           0.62, 1.26, 0.85,
                           1, 4);
    }


    //Spot---------------------------------------
    TEST(runTest_G_pg_0){
        double d = 7.0;
        runSystemTest_G_pg(1, 1.0, 9, 0.1,
                           1.0, d, 0,
                           0.62, 1.26, 0.85,
                           1, 4);
    }


    TEST(runTest_G_pg_1){
        double d = 0.1;
        runSystemTest_G_pg(1, 1.0, 9, 0.1,
                           1.0, d, 0,
                           0.62, 1.26, 0.85,
                           1, 4);
    }

    TEST(runTest_G_pg_2){
        double d = 10;
        runSystemTest_G_pg(1, 1.0, 9, 0.1,
                           1.0, d, 0,
                           0.62, 1.26, 0.85,
                           1, 4);
    }



    TEST(runTest_G_pg_3){
        double d = 0.3;
        runSystemTest_G_pg(3, 0.5, 6, 0.01,
                           1.0, d, 0,
                           0.06, 0.12, 0.8,
                           1, 4);
    }


    TEST(runTest_G_pg_4){
        double d = 1.7;
        runSystemTest_G_pg(3, 0.5, 8, 0.1,
                           0.5, d, 0,
                           0.25, 0.83, 0.8,
                           1, 4);
    }



    TEST(runTest_G_pg_5){
        double d = 10.7;
        runSystemTest_G_pg(1, 1.0, 8, 0.1,
                           0.5, d, 0,
                           0.62, 1.26, 0.85,
                           1, 4);
    }


    TEST(runTest_G_pg_6){
        double d = 0.01;
        runSystemTest_G_pg(1, 1.0, 8, 0.1,
                           1.0, d, 0,
                           0.62, 1.26, 0.85,
                           1, 4);
    }


    TEST(runTest_G_pg_7){
        double d = 7.1;
        runSystemTest_G_pg(1, 1.0, 8, 0.1,
                           1.0, d, 0,
                           0.62, 1.26, 0.85,
                           1, 4);
    }


}



