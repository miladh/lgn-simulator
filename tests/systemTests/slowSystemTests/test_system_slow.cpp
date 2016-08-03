/**********************************************************************
 *  Test: response of ganglion cell, relay cell, interneuron and
 *  cortical cell with patchgrating grating stimulus:
 *
 *               Rg(r,t) = Wg(r,t) * S(r,t)
 *               Ri(r,t) = [Wg(r,t)Kig(r,t) + Kic(r,t)Kcr(r,t)*Wr(r,t)] * S(r,t)
 *               Rc(r,t) = [Wr(r,t)Kcr(r,t)] * S(r,t)
 *               Rr(r,t) = [[W(r,t)Krg(r,t) + Kri(r,t) Wg(r,t)Kig(r,t)]
 *                         /[1 - Kic(r,t)Kcr(r,t)] - Krc(r,t)Kcr(r,t)] ] * S(r,t)
 *
 *  Analytic source: monte carlo integration
 *
 * ********************************************************************/

#include <lgnSimulator.h>
#include <catch.hpp>
#include "test_system_gric_pg_1.h"

using namespace lgnSimulator;

TEST_CASE("system_gric_pg_1 [slow]"){
    /***
    * system: gric
    * stimulus: patch grating
    * test: center relay cell response with
    *       varying mask, spatial frequency,
    *       and feedback weight
    ***/

    string testLabel = "system_gric_pg_1";
    string sourceFilename = "test_system_gric_pg_1";
    test_system_gric_pg_1 test(testLabel,
                                       sourceFilename,
                                       1e4, 5e6, 1e-4);

    test.runTest();
}


TEST_CASE("system_gric_pg_2 [slow]"){
    /***
    * system: gric
    * stimulus: patch grating
    * test: center interneuron response with
    *       varying mask, spatial frequency,
    *       and feedback weight
    ***/

    string testLabel = "system_gric_pg_2";
    string sourceFilename = "test_system_gric_pg_2";
    test_system_gric_pg_1 test(testLabel,
                                       sourceFilename,
                                       1e4, 5e6, 1e-4);

    test.runTest();
}
