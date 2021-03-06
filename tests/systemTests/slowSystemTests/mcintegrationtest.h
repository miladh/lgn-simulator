#ifndef TEST_SYSTEM_MCINTEGRATION_H
#define TEST_SYSTEM_MCINTEGRATION_H

#include <lgnSimulator.h>
#include <catch.hpp>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>


using namespace lgnSimulator;

class MCintegrationTest
{
public:
    MCintegrationTest(string testLabel, string sourceFilename, double preCalls, double calls, double epsilon);

    virtual void runTest() = 0;
    virtual double integrand(double *k, size_t dim, void *params) = 0;
    double computeIntegral(double *xl, double *xu, void *params);
    void display_results (char *title, double result, double error);
    void writeOutputFile();

//    double integrandContainer(double *k, size_t dim, void *params);

    struct Param{
        MCintegrationTest *a;
        const CircleMaskGrating& S;
        const Kernel & Wg;
        const Kernel & Krg;
        const Kernel & Krc;
        const Kernel & Kcr;
        const double peak;
    };

protected:
    string m_outputFilename;
    string m_hashValue;
    int m_ndim=0;

    double m_computed = 0.0;
    double m_peak = 0.0;
    double m_preCalls = 1;
    double m_calls = 1;
    double m_epsilon = 0;
    bool m_computeMC=true;

    vector<double> m_results;
    double m_error = 0.0;




};

#endif // TEST_SYSTEM_MCINTEGRATION_H
