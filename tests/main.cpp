#include <unittest++/UnitTest++.h>
#include <unittest++/Test.h>
#include <unittest++/TestReporterStdout.h>
#include <unittest++/TestRunner.h>

#include "stimuli/patchgrating.h"
#include "spatialKernels/dog.h"
TEST(dog) {


    DOG dog(1.0, 0.25, 0.85, 0.83);
    CHECK_CLOSE(dog.real({0.5, 0.1}), -0.189791527743, 1e-12);
    CHECK_CLOSE(dog.real({1.2, 1.9}), -0.00025733892027, 1e-12);

    CHECK_CLOSE(dog.complex({0.5, 1.1}), 0.316423256919, 1e-12);
    CHECK_CLOSE(dog.complex({1.5, 0.1}), 0.389361200098, 1e-12);

}


//TEST(loopKernel) {
//    Config cfg;
//    cfg.readFile("../../eDOG/tests/configTests.cfg");
//    ImpulseResponse G(&cfg);
//    CHECK_CLOSE(G.loopKernel(0.1, 5.9), 0.00124325581819, 1e-12);
//    CHECK_CLOSE(G.loopKernel(1.5, 0.1), 0.338789712712, 1e-12);

//}

TEST(heaviside) {
    CHECK_EQUAL(Functions::heaviside(-1.2), 0);
    CHECK_EQUAL(Functions::heaviside(2.2), 1.);

}

TEST(stimuli) {
    // TODO: Need complex stimuli tests!

    //    Config cfg;
    //    cfg.readFile("../../eDOG/tests/configTests.cfg");
    //    PatchGrating S(&cfg);
    //    CHECK_CLOSE(S.spatial({-0.1, 0.1}, 0.5), 5.81683089464, 1e-11);
    //    CHECK_CLOSE(S.spatial({-0.7, 0.9}, 2.5), 9.97768728668, 1e-11);


}


int main()
{
    int result = 0;
    bool fft_fft = 0;
    bool dev = 1;


    UnitTest::TestReporterStdout reporter;
    UnitTest::TestRunner runner(reporter);

        if(dev){
      result += runner.RunTestsIf(UnitTest::Test::GetTestList(), "DEVELOPMENT", UnitTest::True(), 0);
        }

    if(fft_fft){
        result += runner.RunTestsIf(UnitTest::Test::GetTestList(), "FFT", UnitTest::True(), 0);
    }
    return result;

    //    return UnitTest::RunAllTests();
}
