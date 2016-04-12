#include <iostream>
#include <unistd.h>
#include <time.h>
#include <yaml-cpp/yaml.h>

#include <lgnSimulator.h>

using namespace std;
using namespace lgnSimulator;

int main(int argc, char* argv[])
{

    cout << "=====LGN Simulator Model: Spatial summation =====" << endl;
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
    io.writeIntegratorProperties(integrator);

    //Stim---------------------------------------------------------------------
    unique_ptr<Grating> S = createGratingStimulus(integrator, cfg["stimulus"]);

    //Ganglion cell:-----------------------------------------------------------
    DOG Wg_s = createSpatialDOGKernel(cfg["ganglion"]["Wg"]);
    Biphasic Wg_t = createTemporalBiphasicKernel(cfg["ganglion"]["Wt"]);

    SeparableKernel Wg(&Wg_s, &Wg_t);
    GanglionCell ganglion(integrator, Wg, cfg["ganglion"]["R0"].as<double>());


    //Relay cell: -------------------------------------------------------------
    RelayCell relay(integrator, cfg["relay"]["R0"].as<double>());

    // G -> R
    SpatialGaussian Ks_rg = createSpatialGaussianKernel(cfg["relay"]["Krg"]["spatial"]);
    TemporalDelta Kt_rg = createTemporalDeltaKernel(cfg["relay"]["Krg"]["temporal"]);
    SeparableKernel Krg(&Ks_rg, &Kt_rg);


    // I -> R
    SpatialGaussian Ks_ri = createSpatialGaussianKernel(cfg["relay"]["Kri"]["spatial"]);
    TemporalDelta Kt_ri = createTemporalDeltaKernel(cfg["relay"]["Kri"]["temporal"]);
    SeparableKernel Kri(&Ks_ri, &Kt_ri);

    // C -> R
    SpatialGaussian Ks_rc = createSpatialGaussianKernel(cfg["relay"]["Krc"]["spatial"]);
    TemporalDelta Kt_rc = createTemporalDeltaKernel(cfg["relay"]["Krc"]["temporal"]);

    //    EllipticGaussian Ks_rc = createSpatialEllipticGaussianKernel(cfg["relay"]["Krc"]["spatial"]));
    //    DOE Kt_rc = createTemporalDOEKernel(cfg["relay"]["Krc"]["temporal"]);

    SeparableKernel Krc(&Ks_rc, &Kt_rc);

    //Relay cell: -------------------------------------------------------------
    Interneuron interneuron(integrator, cfg["interneuron"]["R0"].as<double>());

    // G -> I
    SpatialGaussian Ks_ig = createSpatialGaussianKernel(cfg["interneuron"]["Kig"]["spatial"]);
    TemporalDelta Kt_ig = createTemporalDeltaKernel(cfg["interneuron"]["Kig"]["temporal"]);
    SeparableKernel Kig(&Ks_ig, &Kt_ig);

    // C -> I
    SpatialGaussian Ks_ic = createSpatialGaussianKernel(cfg["interneuron"]["Kic"]["spatial"]);
    TemporalDelta Kt_ic = createTemporalDeltaKernel(cfg["interneuron"]["Kic"]["temporal"]);
    SeparableKernel Kic(&Ks_ic, &Kt_ic);



    //Cortical cell: -------------------------------------------------------------
    CorticalCell cortical(integrator, cfg["cortical"]["R0"].as<double>());

    // R -> G
    SpatialGaussian Ks_cr = createSpatialGaussianKernel(cfg["cortical"]["Kcr"]["spatial"]);
    TemporalDelta Kt_cr = createTemporalDeltaKernel(cfg["cortical"]["Kcr"]["temporal"]);
    //    DOE Kt_cr = createTemporalDOEKernel(cfg["cortical"]["Kcr"]["temporal"]);

    SeparableKernel Kcr(&Ks_cr, &Kt_cr);


    //Connect neurons:---------------------------------------------------------
    relay.addGanglionCell(&ganglion, Krg);
    relay.addCorticalCell(&cortical, Krc);
    relay.addInterNeuron(&interneuron, Kri);

    interneuron.addGanglionCell(&ganglion, Kig);
    interneuron.addCorticalCell(&cortical, Kic);

    cortical.addRelayCell(&relay, Kcr);

    //Compute:-----------------------------------------------------------------

    if(cfg["stimulus"]["storeSpatiotemporal"].as<bool>()){
        S->computeSpatiotemporal();
        io.writeStimulus(S.get(),cfg["stimulus"]["storeFT"].as<bool>());
        S->clearSpatioTemporal();
    }
    S->computeFourierTransform();
    ganglion.computeImpulseResponseFourierTransform();
    relay.computeImpulseResponseFourierTransform();
    cortical.computeImpulseResponseFourierTransform();



    if(cfg["ganglion"]["storeResponse"].as<bool>()){
        ganglion.computeResponse(S.get());
        io.writeResponse(ganglion,
                         cfg["ganglion"]["storeResponseFT"].as<bool>());
        ganglion.clearResponse();
    }

    if( cfg["ganglion"]["storeImpulseResponse"].as<bool>()){
        ganglion.computeImpulseResponse();
        io.writeImpulseResponse(ganglion,
                                cfg["ganglion"]["storeImpulseResponseFT"].as<bool>());
        ganglion.clearImpulseResponse();
    }

    /////////////////////////////////////////////////////////////////////////////////////
    if(cfg["relay"]["storeResponse"].as<bool>()){
        relay.computeResponse(S.get());
        io.writeResponse(relay,
                         cfg["relay"]["storeResponseFT"].as<bool>());
        relay.clearResponse();
    }

    if( cfg["relay"]["storeImpulseResponse"].as<bool>()){
        relay.computeImpulseResponse();
        io.writeImpulseResponse(relay,
                                cfg["relay"]["storeImpulseResponseFT"].as<bool>());
        relay.clearImpulseResponse();
    }


    /////////////////////////////////////////////////////////////////////////////////////
    if(cfg["interneuron"]["storeResponse"].as<bool>()){
        interneuron.computeResponse(S.get());
        io.writeResponse(interneuron,
                         cfg["interneuron"]["storeResponseFT"].as<bool>());
        interneuron.clearResponse();
    }

    if( cfg["interneuron"]["storeImpulseResponse"].as<bool>()){
        interneuron.computeImpulseResponse();
        io.writeImpulseResponse(interneuron,
                                cfg["interneuron"]["storeImpulseResponseFT"].as<bool>());
        interneuron.clearImpulseResponse();
    }

    /////////////////////////////////////////////////////////////////////////////////////
    if(cfg["cortical"]["storeResponse"].as<bool>()){
        cortical.computeResponse(S.get());
        io.writeResponse(cortical,
                         cfg["cortical"]["storeResponseFT"].as<bool>());
        cortical.clearResponse();
    }

    if( cfg["cortical"]["storeImpulseResponse"].as<bool>()){
        cortical.computeImpulseResponse();
        io.writeImpulseResponse(cortical,
                             cfg["cortical"]["storeImpulseResponseFT"].as<bool>());
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

