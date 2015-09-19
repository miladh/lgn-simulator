#include <unittest++/UnitTest++.h>
#include <unittest++/Test.h>
#include <unittest++/TestReporterStdout.h>
#include <unittest++/TestRunner.h>


int main()
{
    int result = 0;

    bool mathTests = 1;
    bool stimTests = 0;
    bool fft_1D = 0;
    bool dev = 1;

    UnitTest::TestReporterStdout reporter;
    UnitTest::TestRunner runner(reporter);

    if(dev){
        result += runner.RunTestsIf(UnitTest::Test::GetTestList(), "DEVELOPMENT", UnitTest::True(), 0);
    }


    if(stimTests){
        result += runner.RunTestsIf(UnitTest::Test::GetTestList(), "STIMULI", UnitTest::True(), 0);
    }

    if(mathTests){
        result += runner.RunTestsIf(UnitTest::Test::GetTestList(), "SPECIALFUNCTIONS", UnitTest::True(), 0);
    }

    if(fft_1D){
        result += runner.RunTestsIf(UnitTest::Test::GetTestList(), "FFT_1D", UnitTest::True(), 0);
    }


    return result;

    //    return UnitTest::RunAllTests();
}
