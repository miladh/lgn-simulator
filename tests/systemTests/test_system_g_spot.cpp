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

double response(double d, double s, double A, double a, double B, double b)
{
    double l = A * s;
    double w = B/A;

    double Rg = l * (  1 - exp(-d*d / (4*a*a)) - w * (1. -  exp(-d*d / (4*b*b) ) ) );

    return Special::heaviside(Rg) * Rg;
}


void runSystemTest_G_spot(int nt, double dt, int ns, double ds,
                          double C, double maskSize, int delayId,
                          double A, double a, double B, double b)
{

    //integrator
    Integrator integrator(nt, dt, ns, ds);
    vec t = integrator.timeVec();

    //stimulus
    CircleMaskGrating spot(integrator, 0, 0, 0, C, maskSize);
    spot.computeFourierTransform();

    //ganglion cell
    DOG Ws(A, a, B, b);
    TemporalDelta Wt(t[delayId], dt);
    SeparableKernel W(&Ws, &Wt);

    GanglionCell ganglion(integrator, W);
    ganglion.computeResponse(&spot);


    //Compute analytic:
    double Rg_ex = response(spot.maskSize(), C,  A, a, B, b);

    //Compute numerical
    ganglion.computeResponse(&spot);
    double Rgc = ganglion.response()(integrator.nPointsSpatial()/2,
                                     integrator.nPointsSpatial()/2,
                                     delayId);


    //Test
    cout << Rg_ex - Rgc << endl;
    CHECK_CLOSE(Rg_ex, Rgc, 1e-9);
}

SUITE(system){
    TEST(runTest_G_spot_0){
        runSystemTest_G_spot(3, 0.5, 9, 0.01,
                            0.1, 0.5, 0,
                            1, 0.06, 0.8, 0.12);
    }
}





