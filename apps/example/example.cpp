#include <iostream>
#include <yaml-cpp/yaml.h>
#include <lgnSimulator.h>

using namespace std;
using namespace lgnSimulator;

int main(int argc, char* argv[])
{
    cout << "===== lgn-simulator =====" << endl;

    //read config file-----------------------------------------------------------------------
    YAML::Node cfg = YAML::LoadFile(argv[1]);

    //Output manager:---------------------------------------------------------
    OutputManager outputManager(cfg["OutputManager"]["outputFilename"].as<std::string>());

    //Integrator-----------------------------------------------------------------------------
    Integrator integrator = createIntegrator(cfg["grid"]);

    //Stim-----------------------------------------------------------------------------------
    unique_ptr<Grating> S = createGratingStimulus(&integrator, cfg["stimulus"]);

    //Ganglion cell:--------------------------------------------------------------------------
    DOG Wg_s = createSpatialDOGKernel(cfg["ganglion"]["Wg"]);
    TemporalDelta Wg_t = createTemporalDeltaKernel(cfg["ganglion"]["Wt"]);
    SeparableKernel Wg(cfg["ganglion"]["w"].as<double>(), &Wg_s, &Wg_t);

    GanglionCell ganglion(&integrator, Wg, cfg["ganglion"]["R0"].as<double>());

    //Relay cell: ------------------------------------------------------------------------------
    RelayCell relay(&integrator, cfg["relay"]["R0"].as<double>());

    // G -> R
    SpatialGaussian Ks_rg = createSpatialGaussianKernel(cfg["relay"]["Krg"]["spatial"]);
    TemporalDelta Kt_rg = createTemporalDeltaKernel(cfg["relay"]["Krg"]["temporal"]);
    SeparableKernel Krg(cfg["relay"]["Krg"]["w"].as<double>(), &Ks_rg, &Kt_rg);

    // C -> R
    SpatialGaussian Ks_rc = createSpatialGaussianKernel(cfg["relay"]["Krc"]["spatial"]);
    TemporalDelta Kt_rc = createTemporalDeltaKernel(cfg["relay"]["Krc"]["temporal"]);
    SeparableKernel Krc(cfg["relay"]["Krc"]["w"].as<double>(), &Ks_rc, &Kt_rc);

    //Cortical cell: -----------------------------------------------------------------------------
    HeavisideNonlinearity heavisideNonlinearity;
    CorticalCell cortical(&integrator, cfg["cortical"]["R0"].as<double>(), &heavisideNonlinearity);

    // R -> C
    SpatialDelta Ks_cr = createSpatialDeltaKernel(cfg["cortical"]["Kcr"]["spatial"]);
    TemporalDelta Kt_cr = createTemporalDeltaKernel(cfg["cortical"]["Kcr"]["temporal"]);
    SeparableKernel Kcr(cfg["cortical"]["Kcr"]["w"].as<double>(), &Ks_cr, &Kt_cr);


    //Connect neurons:------------------------------------------------------------------------------
    relay.addGanglionCell(&ganglion, Krg);
    relay.addCorticalCell(&cortical, Krc);
    cortical.addRelayCell(&relay, Kcr);

    //Compute:--------------------------------------------------------------------------------------
    S->computeFourierTransform();
    relay.computeImpulseResponse();
    relay.computeResponse(S.get());

    //Write:-----------------------------------------------------------------------------------------
    outputManager.writeIntegratorProperties(integrator);
    outputManager.writeStimulusProperties(S.get());
    outputManager.writeStimulus(S.get());
    outputManager.writeResponse(relay);
    outputManager.writeImpulseResponse(relay);

    return 0;
}
