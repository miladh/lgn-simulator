//#include <unittest++/UnitTest++.h>
//#include <unittest++/Test.h>
//#include <unittest++/TestReporterStdout.h>
//#include <unittest++/TestRunner.h>

#define CATCH_CONFIG_MAIN
#include <catch.hpp>
//#include <iostream>


//using namespace UnitTest;
//using namespace std;

//int main()
//{
//    int result = 0;

//    bool special = 1;
//    bool stimulus =1 ;
//    bool kernel = 1;
//    bool fftHelperTests = 1;
//    bool integrator = 1;
//    bool systemTests = 1;

//    TestReporterStdout reporter;
//    TestRunner runner(reporter);

//    if(special){
//        cout << "Running special tests..." << endl;
//        result += runner.RunTestsIf(Test::GetTestList(), "special", True(), 0);
//        cout << "Special tests completed. " << endl << endl;
//    }

//    if(stimulus){
//        cout << "Running stimulus tests..." << endl;
//        result += runner.RunTestsIf(Test::GetTestList(), "stimulus", True(), 0);
//        cout << "Stimulus tests completed. " << endl << endl;
//    }

//    if(kernel){
//        cout << "Running kernel tests..." << endl;
//        result += runner.RunTestsIf(Test::GetTestList(), "kernel", True(), 0);
//        cout << "Kernel tests completed. " << endl << endl;
//    }

//    if(fftHelperTests){
//        cout << "Running FFT helper tests..." << endl;
//        result += runner.RunTestsIf(Test::GetTestList(), "fftHelper", True(), 0);
//        cout << "FFT helper tests completed. " << endl << endl;
//    }

//    if(integrator){
//        cout << "Running integrator tests..." << endl;
//        result += runner.RunTestsIf(Test::GetTestList(), "integrator", True(), 0);
//        cout << "Integrator tests completed. " << endl << endl;
//    }

//    if(systemTests){
//        cout << "Running system tests..." << endl;
//        result += runner.RunTestsIf(Test::GetTestList(), "system", True(), 0);
//        cout << "System tests completed. " << endl << endl;
//    }



//    return result;

//    //    return RunAllTests();
//}
