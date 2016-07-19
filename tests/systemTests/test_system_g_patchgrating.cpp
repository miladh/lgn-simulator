/**********************************************************************
 *  Test: ganglion cell response for patch grating, in the case where
 *        the stimulation spot and receptive-field are concentric.
 *
 *  Analytic source: Einevoll et. al (2005)
 *
 * ********************************************************************/


#include <lgnSimulator.h>
#include <catch.hpp>


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
                        double weight, unsigned int n,
                        double eps)
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
    //    cout << kd << endl;
    //    cout << Rg_ex - Rgc << endl;
    INFO( "kd=" << kd << "  d=" << maskSize);
    REQUIRE(Rg_ex== Approx(Rgc).epsilon(eps));
}





//Grating---------------------------------------------
TEST_CASE("runTest_G_pg_grating"){
    double eps = 1e-3;
    vec diameters = linspace(0.0, 10, 11);
    vec kId =linspace(0, 90, 10);
    for(double d : diameters){
        for(double k : kId){
            runSystemTest_G_pg(1, 1.0, 9, 0.1,
                               1.0, d, k,
                               0.62, 1.26, 0.85,
                               1, 4, eps);
        }
    }

}


//Spot--------------------------------------------
TEST_CASE("runTest_G_pg_spot"){
    double eps = 1e-6;
    vec diameters = linspace(0.0, 10, 11);
    for(double d : diameters){
            runSystemTest_G_pg(1, 1.0, 9, 0.1,
                               1.0, d, 0,
                               0.62, 1.26, 0.85,
                               1, 4, eps);
    }



}





