#include <unittest++/UnitTest++.h>
#include <unittest++/Test.h>
#include <unittest++/TestReporterStdout.h>
#include <unittest++/TestRunner.h>



int main()
{
    int result = 0;

    bool mathTests = 1;
    bool stimTests = 1;
    bool fftHelperTests = 1;
    bool integratorTests = 1;
    bool systemTests = 1;

    bool dev = 1;


    UnitTest::TestReporterStdout reporter;
    UnitTest::TestRunner runner(reporter);

    if(stimTests){
        result += runner.RunTestsIf(UnitTest::Test::GetTestList(), "STIMULI", UnitTest::True(), 0);
    }

    if(mathTests){
        result += runner.RunTestsIf(UnitTest::Test::GetTestList(), "SPECIALFUNCTIONS", UnitTest::True(), 0);
    }

    if(fftHelperTests){
        result += runner.RunTestsIf(UnitTest::Test::GetTestList(), "fftHelper", UnitTest::True(), 0);
    }

    if(integratorTests){
        result += runner.RunTestsIf(UnitTest::Test::GetTestList(), "INTEGRATOR", UnitTest::True(), 0);
    }



    if(dev){
        result += runner.RunTestsIf(UnitTest::Test::GetTestList(), "DEVELOPMENT", UnitTest::True(), 0);
    }


    if(systemTests){
        result += runner.RunTestsIf(UnitTest::Test::GetTestList(), "SYSTEM", UnitTest::True(), 0);
    }



    return result;

    //    return UnitTest::RunAllTests();
}
