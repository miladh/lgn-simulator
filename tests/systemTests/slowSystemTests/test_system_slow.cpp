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
#include "test_system_gric_patchgrating.h"

using namespace lgnSimulator;

TEST_CASE("system_gric_pg_1 [slow]"){
    /***
    * system: gric
    * stimulus: patch grating
    * test: center relay cell response with
    *       varying mask and spatial frequency
    ***/

    test_system_girc_patchGrating test("system_gric_pg",
                                       "test_system_gric_patchgrating",
                                       1e3, 5e6);

    test.runTest();
}
