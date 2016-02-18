/**********************************************************************
 *  Test: various mathematical helper functions
 *
 *  Analytic source: by hand and Python
 *
 * ********************************************************************/

#include <unittest++/UnitTest++.h>
#include "kernels/spatialKernels/dog.h"


using namespace lgnSimulator;

SUITE(SPECIALFUNCTIONS){


    TEST(secondKindBesselFunction) {
        CHECK_EQUAL(Functions::secondKindBesselFunction(0.0), 0);
        CHECK_CLOSE(Functions::secondKindBesselFunction(0.1), 0.049937526036241998, 1e-12);
        CHECK_CLOSE(Functions::secondKindBesselFunction(5.0), -0.32757913759146522, 1e-12);
        CHECK_CLOSE(Functions::secondKindBesselFunction(7.01558666981562), 0.0, 1e-12);
        CHECK_CLOSE(Functions::secondKindBesselFunction(-7.01558666981562), 0.0, 1e-12);
        CHECK_CLOSE(Functions::secondKindBesselFunction(3.83170597020751), 0.0, 1e-12);
    }


    TEST(heaviside) {
        CHECK_EQUAL(Functions::heaviside(-1.2), 0);
        CHECK_EQUAL(Functions::heaviside(0.0), 1.);
        CHECK_EQUAL(Functions::heaviside(2.2), 1.);
    }

    TEST(delta) {
        CHECK_EQUAL(Functions::delta(0,0), 1);
        CHECK_EQUAL(Functions::delta(-1e6,-1e6), 1);
        CHECK_EQUAL(Functions::delta(3.456,3.456), 1);
        CHECK_EQUAL(Functions::delta(2.1,-2.1), 0);
        CHECK_EQUAL(Functions::delta(1.,3.), 0);
        CHECK_EQUAL(Functions::delta(-1230.,-50.), 0);
    }

    TEST(isOdd) {
        CHECK_EQUAL(Functions::isOdd(0), 0);
        CHECK_EQUAL(Functions::isOdd(-3),1);
        CHECK_EQUAL(Functions::isOdd(2), 0);
        CHECK_EQUAL(Functions::isOdd(2.9), 0);
        CHECK_EQUAL(Functions::isOdd(-3.6),1);
        CHECK_EQUAL(Functions::isOdd(91), 1);
    }


    TEST(dog) {
        DOG dog(1.0, 0.25, 0.85, 0.83);
        CHECK_CLOSE(dog.spatial({0.5, 0.1}), -0.189791527743, 1e-12);
        CHECK_CLOSE(dog.spatial({1.2, 1.9}), -0.00025733892027, 1e-12);

        CHECK_CLOSE(real(dog.fourierTransform({0.5, 1.1})), 0.316423256919, 1e-12);
        CHECK_CLOSE(real(dog.fourierTransform({1.5, 0.1})), 0.389361200098, 1e-12);

        CHECK_EQUAL(imag(dog.fourierTransform({0.5, 1.1})), 0.0);
        CHECK_EQUAL(imag(dog.fourierTransform({1.5, 0.1})), 0.0);

    }




}
