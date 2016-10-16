#include <iostream>
#include <unistd.h>
#include <time.h>
#include <yaml-cpp/yaml.h>

#include <lgnSimulator.h>

using namespace std;
using namespace lgnSimulator;

int main(int argc, char* argv[])
{

    cout << "===== LGN Simulator Model: spatiotemporal extended DoG =====" << endl;
    clock_t t;
    t = clock();

    if(argc < 2) {
        cerr << "Too few arguments." << endl;
        return 1;
    }


    //read config file-------------------------------------------------------
    YAML::Node cfg = YAML::LoadFile(argv[1]);

    //Output manager:---------------------------------------------------------
    OutputManager io(cfg["OutputManager"]["outputFilename"].as<std::string>());

    //Integrator--------------------------------------------------------------
    Integrator integrator = createIntegrator(cfg["grid"]);

    //Stim---------------------------------------------------------------------
    unique_ptr<Grating> S = createGratingStimulus(&integrator, cfg["stimulus"]);

    //Ganglion cell:-----------------------------------------------------------
    DOG Wg_s = createSpatialDOGKernel(cfg["ganglion"]["Wg"]);
    Biphasic Wg_t = createTemporalBiphasicKernel(cfg["ganglion"]["Wt"]);
    SeparableKernel Wg(cfg["ganglion"]["w"].as<double>(), &Wg_s, &Wg_t);

    GanglionCell ganglion(&integrator, Wg, cfg["ganglion"]["R0"].as<double>());


    //Relay cell: -------------------------------------------------------------
    RelayCell relay(&integrator, cfg["relay"]["R0"].as<double>());

    // G -> R
    SpatialGaussian Ks_rg = createSpatialGaussianKernel(cfg["relay"]["Krg"]["spatial"]);
    DecayingExponential Kt_rg = createTemporalDecayingExponentialKernel(cfg["relay"]["Krg"]["temporal"]);
    SeparableKernel Krg(cfg["relay"]["Krg"]["w"].as<double>(), &Ks_rg, &Kt_rg);

    // G -> I -> R
    SpatialGaussian Ks_rig = createSpatialGaussianKernel(cfg["relay"]["Krig"]["spatial"]);
    DecayingExponential Kt_rig = createTemporalDecayingExponentialKernel(cfg["relay"]["Krig"]["temporal"]);
    SeparableKernel Krig(cfg["relay"]["Krig"]["w"].as<double>(), &Ks_rig, &Kt_rig);

    // C -> R
    DOG Ks_rc = createSpatialDOGKernel(cfg["relay"]["Krc"]["spatial"]);
    DecayingExponential Kt_rc = createTemporalDecayingExponentialKernel(cfg["relay"]["Krc"]["temporal"]);
    SeparableKernel Krc(cfg["relay"]["Krc"]["w"].as<double>(), &Ks_rc, &Kt_rc);


    //Cortical cell: -------------------------------------------------------------
    HeavisideNonlinearity heavisideNonlinearity;
    CorticalCell cortical(&integrator, cfg["cortical"]["R0"].as<double>(), &heavisideNonlinearity);

    // R -> C
    SpatialDelta Ks_cr = createSpatialDeltaKernel(cfg["cortical"]["Kcr"]["spatial"]);
    TemporalDelta Kt_cr = createTemporalDeltaKernel(cfg["cortical"]["Kcr"]["temporal"]);
    SeparableKernel Kcr(cfg["cortical"]["Kcr"]["w"].as<double>(), &Ks_cr, &Kt_cr);


    //Connect neurons:---------------------------------------------------------
    relay.addGanglionCell(&ganglion, Krg);
    relay.addGanglionCell(&ganglion, Krig);
    relay.addCorticalCell(&cortical, Krc);
    cortical.addRelayCell(&relay, Kcr);

    //Compute:-----------------------------------------------------------------
    S->computeFourierTransform();
    ganglion.computeImpulseResponseFourierTransform();
    relay.computeImpulseResponseFourierTransform();

    //Write:-----------------------------------------------------------------
    io.copyConfigFile(argv[1]);
    io.writeIntegratorProperties(integrator);
    io.writeStimulusProperties(S.get());

    /////////////////////////////////////////////////////////////////////////////////////
    if(cfg["stimulus"]["storeFT"].as<bool>()){
        io.writeStimulusFourierTransform(S.get());
    }
    if(cfg["stimulus"]["storeSpatiotemporal"].as<bool>()){
        S->computeSpatiotemporal();
        io.writeStimulus(S.get());
        S->clearSpatioTemporal();
    }

    /////////////////////////////////////////////////////////////////////////////////////
    if(cfg["ganglion"]["storeResponseFT"].as<bool>()){
        ganglion.computeResponseFourierTransform(S.get());
        io.writeResponseFourierTransform(ganglion);
    }
    if(cfg["ganglion"]["storeResponse"].as<bool>()){
        ganglion.computeResponse(S.get());
        io.writeResponse(ganglion);
        ganglion.clearResponse();
    }
    ganglion.clearResponseFourierTransform();
    if(cfg["ganglion"]["storeImpulseResponseFT"].as<bool>()){
        io.writeImpulseResponseFourierTransform(ganglion);
    }
    if( cfg["ganglion"]["storeImpulseResponse"].as<bool>()){
        ganglion.computeImpulseResponse();
        io.writeImpulseResponse(ganglion);
        ganglion.clearImpulseResponse();
    }

    /////////////////////////////////////////////////////////////////////////////////////
    if(cfg["relay"]["storeResponseFT"].as<bool>()){
        relay.computeResponseFourierTransform(S.get());
        io.writeResponseFourierTransform(relay);
    }
    if(cfg["relay"]["storeResponse"].as<bool>()){
        relay.computeResponse(S.get());
        io.writeResponse(relay);
        relay.clearResponse();
    }
    relay.clearResponseFourierTransform();

    if( cfg["relay"]["storeImpulseResponseFT"].as<bool>()){
        io.writeImpulseResponseFourierTransform(relay);
    }
    if( cfg["relay"]["storeImpulseResponse"].as<bool>()){
        relay.computeImpulseResponse();
        io.writeImpulseResponse(relay);
        relay.clearImpulseResponse();
    }


    /////////////////////////////////////////////////////////////////////////////////////
    if(cfg["cortical"]["storeResponseFT"].as<bool>()){
        cortical.computeResponseFourierTransform(S.get());
        io.writeResponseFourierTransform(cortical);
    }
    if(cfg["cortical"]["storeResponse"].as<bool>()){
        cortical.computeResponse(S.get());
        io.writeResponse(cortical);
        cortical.clearResponse();
    }
    cortical.clearResponseFourierTransform();

    if( cfg["cortical"]["storeImpulseResponseFT"].as<bool>()){
        io.writeImpulseResponseFourierTransform(cortical);
    }
    if( cfg["cortical"]["storeImpulseResponse"].as<bool>()){
        cortical.computeImpulseResponse();
        io.writeImpulseResponse(cortical);
        cortical.clearImpulseResponse();
    }


    //Finalize:----------------------------------------------------------
    t = clock() - t;
    double elapsedTime = ((float)t)/CLOCKS_PER_SEC;
    if(elapsedTime <= 60){
        printf ("%f seconds.\n", elapsedTime);
    }else{
        printf ("%f minutes.\n", elapsedTime/60);
    }


    return 0;
}
