#ifndef TEST_SYSTEM_GIRC_PG_1_H
#define TEST_SYSTEM_GIRC_PG_1_H

#include "mcintegrationtest.h"


class test_system_grc_pg_1 : public MCintegrationTest
{
public:
    test_system_grc_pg_1(string testLabel, string filename,
                                  double preCalls, double calls, double epsilon);

    // test_system_MCintegration interface
public:
    virtual void runTest() override;
    virtual double integrand(double *k, size_t dim, void *params) override;
};

#endif // TEST_SYSTEM_GIRC_PG_1_H
