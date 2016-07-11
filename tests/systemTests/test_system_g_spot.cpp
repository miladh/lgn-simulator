/**********************************************************************
 *  Test: ganglion cell response for a static spot, in the case where
 *        the stimulation spot and receptive-field are concentric.
 *
 *  Analytic source: Einevoll et. al (2000)
 *
 * ********************************************************************/


#include <unittest++/UnitTest++.h>
#include <lgnSimulator.h>

using namespace lgnSimulator;

double response(double d, double s, double a, double b, double c)
{

    double Rg = s * (  1 - exp(-d*d / (4*a*a)) - c * (1. -  exp(-d*d / (4*b*b) ) ) );

    return Special::heaviside(Rg) * Rg;
}


void runSystemTest_G_spot(int nt, double dt, int ns, double ds,
                          double C, double maskSize, int delayId,
                          double a, double b, double c,
                          double weight)
{

    //integrator
    Integrator integrator(nt, dt, ns, ds);
    vec t = integrator.timeVec();

    //stimulus
    CircleMaskGrating spot(integrator, 0, 0, 0, C, maskSize);
    spot.computeFourierTransform();

    //ganglion cell
    DOG Ws(a, b, c);
    TemporalDelta Wt(t[delayId], dt);
    SeparableKernel W(weight, &Ws, &Wt);

    GanglionCell ganglion(integrator, W);
    ganglion.computeResponse(&spot);


    //Compute analytic:
    double Rg_ex = response(spot.maskSize(), C, a, b, c);

    //Compute numerical
    ganglion.computeResponse(&spot);
    double Rgc = ganglion.response()(integrator.nPointsSpatial()/2,
                                     integrator.nPointsSpatial()/2,
                                     delayId);


    //Test
//    cout << Rg_ex - Rgc << endl;
    CHECK_CLOSE(Rg_ex, Rgc, 1e-8);
}

SUITE(system){
    TEST(runTest_G_spot_0){
        runSystemTest_G_spot(3, 0.5, 6, 0.01,
                            0.1, 0.3, 0,
                            0.06, 0.12, 0.8,
                             1);
    }

    TEST(runTest_G_spot_1){
        runSystemTest_G_spot(3, 0.5, 8, 0.1,
                            0.5, 1.7, 4,
                            0.25, 0.83, 0.8,
                            1);
    }

    TEST(runTest_G_spot_2){
        runSystemTest_G_spot(1, 1.0, 8, 0.1,
                            1.0, 14.7, 0,
                            0.62, 1.26, 0.85,
                            1);
    }

    TEST(runTest_G_spot_3){
        runSystemTest_G_spot(1, 1.0, 8, 0.1,
                            1.0, 0.01, 0,
                            0.62, 1.26, 0.85,
                            1);
    }
}





