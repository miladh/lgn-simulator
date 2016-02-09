/**********************************************************************
 *  Test: ganglion cell response for a static spot, in the case where
 *        the stimulation spot and receptive-field are concentric.
 *
 *  Analytic source: Einevoll et. al (2000)
 *
 * ********************************************************************/

#include <unittest++/UnitTest++.h>
#include "integrator.h"
#include "neurons/ganglioncell.h"
#include "spatialKernels/dog.h"
#include "temporalKernels/temporallyconstant.h"
#include "stimuli/grating/circlemaskgrating.h"

using namespace lgnSimulator;

double response(double d, double s, double A, double a, double B, double b)
{
    double l = A * s;
    double w = B/A;

    double Rg = l * (  1 - exp(-d*d / (4*a*a)) - w * (1. -  exp(-d*d / (4*b*b) ) ) );

    return Rg;
}

SUITE(SYSTEM){


    TEST(ganglion_response_spot){
        int ns = 9;
        int nt = 3;
        double dt = 1.2;
        double ds = 0.01;

        double C = 2.0;
        double maskSize = 0.9;

        double A = 1.0;
        double a = 0.06;
        double B = 0.8;
        double b = 0.12;
        double t0 = 3.0;

        int Ns = pow(2,ns);

        //Integrator
        Integrator integrator(nt, dt, ns, ds);

        //Stimulus
        CircleMaskGrating S(&integrator, {0, 0}, 0, C, maskSize);
        S.computeFourierTransform();

        //Kernels
        DOG Ks(A, a, B, b);
        TemporallyConstant Kt(t0);

        //Cell
        GanglionCell ganglion(&integrator, &Ks, &Kt);

        //Compute analytic:
        double Rg_ex = response(S.maskSize(), t0*C,  A, a, B, b) ;

        //Compute numerical
        ganglion.computeResponse(&S);
        double Rgc = real(ganglion.response()(Ns/2, Ns/2, 5));

        //Test
        CHECK_CLOSE(Rg_ex, Rgc, 1e-9);

    }


}





