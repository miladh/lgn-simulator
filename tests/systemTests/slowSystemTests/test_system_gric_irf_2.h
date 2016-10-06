#ifndef TEST_SYSTEM_GRIC_IRF_2_H
#define TEST_SYSTEM_GRIC_IRF_2_H

#include "mcintegrationtest.h"



class test_system_gric_irf_2 : public MCintegrationTest
{
public:
    test_system_gric_irf_2(string testLabel, string filename, double preCalls, double calls, double epsilon);

    // MCintegrationTest interface
public:
    virtual void runTest() override;
    virtual double integrand(double *k, size_t dim, void *params) override;

private:
    double m_t=0;
};

#endif // TEST_SYSTEM_GRIC_IRF_2_H
