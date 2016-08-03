#ifndef TEST_SYSTEM_GIRC_PATCHGRATING_H
#define TEST_SYSTEM_GIRC_PATCHGRATING_H

#include "mcintegrationtest.h"


class test_system_girc_patchGrating : public MCintegrationTest
{
public:
    test_system_girc_patchGrating(string testLabel, string filename,
                                  double preCalls, double calls);

    // test_system_MCintegration interface
public:
    virtual void runTest() override;
    virtual double integrand(double *k, size_t dim, void *params) override;
};

#endif // TEST_SYSTEM_GIRC_PATCHGRATING_H
