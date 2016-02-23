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
