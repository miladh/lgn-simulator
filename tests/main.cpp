#include <unittest++/UnitTest++.h>
#include <impulseResponse.h>
#include "stimuli/patchgrating.h"
#include "math/dog.h"
#include <lib.h>

TEST(dog) {

    Config cfg;
    cfg.readFile("../../eDOG/tests/configTests.cfg");
    ImpulseResponse G(&cfg);
    CHECK_CLOSE(G.differenceOfGaussian(0.5, 0.1), -0.189791527743, 1e-12);
    CHECK_CLOSE(G.differenceOfGaussian(1.2, 1.9), -0.00025733892027, 1e-12);

    CHECK_CLOSE(G.differenceOfGaussianComplex(0.5, 1.1), 0.316423256919, 1e-12);
    CHECK_CLOSE(G.differenceOfGaussianComplex(1.5, 0.1), 0.389361200098, 1e-12);


    DOG dog(1.0, 0.25, 0.85, 0.83);
    CHECK_CLOSE(dog.real({0.5, 0.1}), -0.189791527743, 1e-12);
    CHECK_CLOSE(dog.real({1.2, 1.9}), -0.00025733892027, 1e-12);

    CHECK_CLOSE(dog.complex({0.5, 1.1}), 0.316423256919, 1e-12);
    CHECK_CLOSE(dog.complex({1.5, 0.1}), 0.389361200098, 1e-12);

}


TEST(loopKernel) {
    Config cfg;
    cfg.readFile("../../eDOG/tests/configTests.cfg");
    ImpulseResponse G(&cfg);
    CHECK_CLOSE(G.loopKernel(0.1, 5.9), 0.00124325581819, 1e-12);
    CHECK_CLOSE(G.loopKernel(1.5, 0.1), 0.338789712712, 1e-12);

}

TEST(heaviside) {
    Functions func;
    CHECK_EQUAL(func.heaviside(-1.2), 0);
    CHECK_EQUAL(func.heaviside(2.2), 1.);

}

TEST(stimuli) {
    Config cfg;
    cfg.readFile("../../eDOG/tests/configTests.cfg");
    PatchGrating S(&cfg);
    CHECK_CLOSE(S.real({-0.1, 0.1}, 0.5), 5.81683089464, 1e-11);
    CHECK_CLOSE(S.real({-0.7, 0.9}, 2.5), 9.97768728668, 1e-11);


}


int main()
{
    return UnitTest::RunAllTests();
}
