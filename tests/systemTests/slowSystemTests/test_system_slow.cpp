/**********************************************************************
 *  Test: response/irf of ganglion cell, relay cell, and
 *  cortical cell with patchgrating grating stimulus:
 *
 *               Rg(r,t) = Wg(r,t) * S(r,t)
 *               Rc(r,t) = [Wr(r,t)Kcr(r,t)] * S(r,t)
 *               Rr(r,t) = [[W(r,t)Krg(r,t)]/[1 - Krc(r,t)Kcr(r,t)] ] * S(r,t)
 *
 *  Analytic source: monte carlo integration
 *
 * ********************************************************************/

#include <lgnSimulator.h>
#include <catch.hpp>
#include "test_system_grc_pg_1.h"
#include "test_system_grc_irf_1.h"
#include "test_system_grc_irf_2.h"

using namespace lgnSimulator;

TEST_CASE("system_grc_pg_1", "[slow]"){
    /***
    * system: grc
    * stimulus: patch grating
    * test: center relay cell response with
    *       varying mask, spatial frequency,
    *       and feedback weight
    ***/

    string testLabel = "system_grc_pg_1";
    string sourceFilename = "test_system_grc_pg_1";
    test_system_grc_pg_1 test(testLabel,
                                       sourceFilename,
                                       1e4, 5e6, 1e-4);

    test.runTest();
}




TEST_CASE("system_grc_irf_1", "[slow]"){
    /***
    * system: grc
    * stimulus: unused
    * test: center relay irf with
    *       varying connectivity weights and size
    ***/

    string testLabel = "system_grc_irf_1";
    string sourceFilename = "test_system_grc_irf_1";
    test_system_grc_irf_1 test(testLabel,
                                       sourceFilename,
                                       1e4, 5e6, 1e-4);

    test.runTest();
}

TEST_CASE("system_grc_irf_2", "[slow]"){
    /***
    * system: grc
    * stimulus: unused
    * test: center relay irf with
    *       including temporal kernels
    ***/

    string testLabel = "system_grc_irf_2";
    string sourceFilename = "test_system_grc_irf_2";
    test_system_grc_irf_2 test(testLabel,
                                       sourceFilename,
                                       1e4, 5e6, 1e-3);

    test.runTest();
}
