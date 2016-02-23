/**********************************************************************
 *  Test: spcial mathematical functions
 *
 *  Analytic source: by hand and Python
 *
 * ********************************************************************/

#include <unittest++/UnitTest++.h>
#include <lgnSimulator.h>

using namespace lgnSimulator;

SUITE(special){


    TEST(secondKindBesselFunction) {
        CHECK_EQUAL(Special::secondKindBesselFunction(0.0), 0);
        CHECK_CLOSE(Special::secondKindBesselFunction(0.1), 0.049937526036241998, 1e-12);
        CHECK_CLOSE(Special::secondKindBesselFunction(5.0), -0.32757913759146522, 1e-12);
        CHECK_CLOSE(Special::secondKindBesselFunction(7.01558666981562), 0.0, 1e-12);
        CHECK_CLOSE(Special::secondKindBesselFunction(-7.01558666981562), 0.0, 1e-12);
        CHECK_CLOSE(Special::secondKindBesselFunction(3.83170597020751), 0.0, 1e-12);
    }


    TEST(heaviside) {
        CHECK_EQUAL(Special::heaviside(-1.2), 0);
        CHECK_EQUAL(Special::heaviside(0.0), 1.);
        CHECK_EQUAL(Special::heaviside(2.2), 1.);
        CHECK_EQUAL(Special::heaviside(-100), 0);
        CHECK_EQUAL(Special::heaviside(1e-6), 1.);
        CHECK_EQUAL(Special::heaviside(34), 1.);
    }

    TEST(delta) {
        CHECK_EQUAL(Special::delta(0,0), 1);
        CHECK_EQUAL(Special::delta(-1e6,-1e6), 1);
        CHECK_EQUAL(Special::delta(3.456,3.456), 1);
        CHECK_EQUAL(Special::delta(2.1,-2.1), 0);
        CHECK_EQUAL(Special::delta(1.,3.), 0);
        CHECK_EQUAL(Special::delta(-1230.,-50.), 0);
    }

    TEST(isOdd) {
        CHECK_EQUAL(Special::isOdd(0), 0);
        CHECK_EQUAL(Special::isOdd(-3),1);
        CHECK_EQUAL(Special::isOdd(2), 0);
        CHECK_EQUAL(Special::isOdd(2.9), 0);
        CHECK_EQUAL(Special::isOdd(-3.6),1);
        CHECK_EQUAL(Special::isOdd(91), 1);
    }



}
