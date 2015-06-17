#include <unittest++/UnitTest++.h>
#include <impulseResponse.h>
#include <stimuli.h>

TEST(dog) {

    Config cfg;
    cfg.readFile("../../eDOG/tests/configTests.cfg");
    ImpulseResponse G(&cfg);
    CHECK_CLOSE(G.differenceOfGaussian(0.5, 0.1), -0.189791527743, 1e-12);
    CHECK_CLOSE(G.differenceOfGaussian(1.2, 1.9), -0.00025733892027, 1e-12);

    CHECK_CLOSE(G.differenceOfGaussianFT(0.5, 1.1), 0.316423256919, 1e-12);
    CHECK_CLOSE(G.differenceOfGaussianFT(1.5, 0.1), 0.389361200098, 1e-12);


}


TEST(loopKernel) {
    Config cfg;
    cfg.readFile("../../eDOG/tests/configTests.cfg");
    ImpulseResponse G(&cfg);
    CHECK_CLOSE(G.loopKernel(0.1, 5.9), 0.00124325581819, 1e-12);
    CHECK_CLOSE(G.loopKernel(1.5, 0.1), 0.338789712712, 1e-12);

}

TEST(heaviside) {
    Config cfg;
    cfg.readFile("../../eDOG/tests/configTests.cfg");
    Stimuli S(&cfg);
    CHECK_EQUAL(S.heaviside(-1.2), 0);
    CHECK_EQUAL(S.heaviside(2.2), 1.);

}

TEST(stimuli) {
    Config cfg;
    cfg.readFile("../../eDOG/tests/configTests.cfg");
    Stimuli S(&cfg);
    CHECK_CLOSE(S.patchGrating(-0.1, 0.1, 0.5), 5.81683089464, 1e-11);
    CHECK_CLOSE(S.patchGrating(-0.7, 0.9, 2.5), 9.97768728668, 1e-11);


}

int main()
{
    return UnitTest::RunAllTests();
}
