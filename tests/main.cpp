#include <unittest++/UnitTest++.h>
#include <unittest++/Test.h>
#include <unittest++/TestReporterStdout.h>
#include <unittest++/TestRunner.h>
#include <iostream>

using namespace UnitTest;
using namespace std;

int main()
{
    int result = 0;


    bool special = 1;
    bool stimulus = 1;
    bool fftHelperTests = 1;
    bool integratorTests = 1;
    bool systemTests = 1;

    bool dev = 1;

    TestReporterStdout reporter;
    TestRunner runner(reporter);

    if(special){
        cout << "Running special tests..." << endl;
        result += runner.RunTestsIf(Test::GetTestList(), "special", True(), 0);
        cout << "Special tests completed. " << endl << endl;
    }

    if(stimulus){
        cout << "Running stimulus tests..." << endl;
        result += runner.RunTestsIf(Test::GetTestList(), "stimulus", True(), 0);
        cout << "Stimulus tests completed. " << endl << endl;
    }


    if(fftHelperTests){
        result += runner.RunTestsIf(Test::GetTestList(), "fftHelper", True(), 0);
    }

    if(integratorTests){
        result += runner.RunTestsIf(Test::GetTestList(), "INTEGRATOR", True(), 0);
    }



    if(dev){
        result += runner.RunTestsIf(Test::GetTestList(), "DEVELOPMENT", True(), 0);
    }


    if(systemTests){
        result += runner.RunTestsIf(Test::GetTestList(), "SYSTEM", True(), 0);
    }



    return result;

    //    return RunAllTests();
}
