/**********************************************************************
 *  Test: spcial mathematical functions
 *
 *  Analytic source: by hand and Python
 *
 * ********************************************************************/

#include <unittest++/UnitTest++.h>
#include <lgnSimulator.h>

using namespace lgnSimulator;

SUITE(stimulus){
    TEST(stimuli) {
        int nt = 5;
        int ns = 5;
        double dt = 0.1;
        double ds = 0.1;
        Integrator integrator(nt, dt, ns, ds);

    }

}
