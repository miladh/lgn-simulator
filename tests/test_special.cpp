/**********************************************************************
 *  Test: spcial mathematical functions
 *
 *  Analytic source: by hand
 *
 * ********************************************************************/

#include <unittest++/UnitTest++.h>
#include <lgnSimulator.h>

using namespace lgnSimulator;

SUITE(special){

    TEST(nearestValue) {
          vec A = {0,1,2,3,4,5};
          double v_A = 2.3;
          CHECK_EQUAL(Special::nearestValue(A, v_A), 2);

          vec B = {-2.3, 5.2, 8.2, 1.1, 2.67,-3.2};
          double v_B = -2.4;
          CHECK_EQUAL(Special::nearestValue(B, v_B), -2.3);

          vec C = {9.66407446,  7.40204369,  8.47934683,  8.7268378 ,  9.537069  ,
                   8.94007828,  6.37876932,  7.84503963,  8.70901142};
          double v_C = 7.5;
          CHECK_EQUAL(Special::nearestValue(C, v_C), 7.40204369);


          vec D = { 4.20844744,  5.44088512, -1.44998235,  1.8764609 , -2.22633141,
                    0.33623971,  7.23507673};
          double v_D = 0.0;
          CHECK_EQUAL(Special::nearestValue(D, v_D), 0.33623971);
    }


    TEST(secondKindBesselFunction) {
        CHECK_EQUAL(Special::secondKindBessel(0.0), 0);
        CHECK_CLOSE(Special::secondKindBessel(0.1), 0.049937526036241998, 1e-12);
        CHECK_CLOSE(Special::secondKindBessel(5.0), -0.32757913759146522, 1e-12);
        CHECK_CLOSE(Special::secondKindBessel(7.01558666981562), 0.0, 1e-12);
        CHECK_CLOSE(Special::secondKindBessel(-7.01558666981562), 0.0, 1e-12);
        CHECK_CLOSE(Special::secondKindBessel(3.83170597020751), 0.0, 1e-12);
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

    TEST(deltaVec2) {
        CHECK_EQUAL(Special::delta({0,0},{0,0}), 1);
        CHECK_EQUAL(Special::delta({-1e6, -1e6}, {-1e6, -1e6}), 1);
        CHECK_EQUAL(Special::delta({3.456,3.456}, {3.456,3.456}), 1);
        CHECK_EQUAL(Special::delta({2.1, 10.78},{2.1, 10.78}), 1);
        CHECK_EQUAL(Special::delta({10.78, 2.1},{2.1, 10.78}), 0);
        CHECK_EQUAL(Special::delta({-4.5, 900.89},{0.2, 67.6}), 0);
        CHECK_EQUAL(Special::delta({2, -2},{2.0001, 2}), 0);
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
