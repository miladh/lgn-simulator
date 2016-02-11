#include <iostream>
#include <unistd.h>
#include <time.h>
#include <yaml-cpp/yaml.h>

#include <lgnSimulator.h>

using namespace std;
using namespace lgnSimulator;

int main(int argc, char* argv[]){

    cout << "========= LGN Simulator: Firing Synchrony =========" << endl;
    clock_t t;
    t = clock();

    if(argc < 2) {
        cerr << "Too few arguments." << endl;
        return 1;
    }

    //read config file-------------------------------------------------------
    YAML::Node cfg = YAML::LoadFile(argv[1]);
    const YAML::Node& ganglionImpRes = cfg["ganglionImpRes"];
    const YAML::Node& Ks_rgSettings = cfg["spatialKernels"]["Krg"];
    const YAML::Node& Ks_rcSettings = cfg["spatialKernels"]["Krc"];
    const YAML::Node& Ks_crSettings = cfg["spatialKernels"]["Kcr"];

    const YAML::Node& Kt_rgSettings = cfg["temporalKernels"]["Krg"];
    const YAML::Node& Kt_rcSettings = cfg["temporalKernels"]["Krc"];
    const YAML::Node& Kt_crSettings = cfg["temporalKernels"]["Kcr"];

    string outputFilename = cfg["outputFile"].as<std::string>();

    //Output manager:---------------------------------------------------------
    OutputManager io(&outputFilename);

    //Integrator--------------------------------------------------------------
    Integrator integrator = createIntegrator(&cfg);

    //Stim---------------------------------------------------------------------
    unique_ptr<Grating> S = createGratingStimulus(&integrator, &cfg);

    //Ganglion cell:-----------------------------------------------------------
    DOG Wg_s = createDOGSpatialKernel(&ganglionImpRes);
    TemporalDelta Wg_t = createTemporalDeltaKernel(&ganglionImpRes);
    GanglionCell ganglion(&integrator, &Wg_s, &Wg_t);

    //Relay cell: -------------------------------------------------------------
    RelayCell relay(&integrator);
    CorticalCell cortical(&integrator);

    //Spatial kernels:---------------------------------------------------------
    SpatialDelta Ks_rg = createSpatialDeltaKernel(&Ks_rgSettings);
    SpatialDelta Ks_cr = createSpatialDeltaKernel(&Ks_crSettings);

    SpatialDelta Ks_rc = createSpatialDeltaKernel(&Ks_rcSettings);
//    Gaussian Ks_rc = createGaussianSpatialKernel(&Ks_rcSettings);

    //Temporal kernels:--------------------------------------------------------
    TemporalDelta Kt_rg = createTemporalDeltaKernel(&Kt_rgSettings);
    TemporalDelta Kt_cr = createTemporalDeltaKernel(&Kt_crSettings);

//    TemporalDelta Kt_rc = createTemporalDeltaKernel(&Kt_rcSettings);
    TemporalGaussian Kt_rc = createGaussianTemporalKernel(&Kt_rcSettings);

    //Connect neurons:---------------------------------------------------------
    relay.addGanglionCell(&ganglion, &Ks_rg, &Kt_rg);
    relay.addCorticalNeuron(&cortical, &Ks_rc, &Kt_rc);
    cortical.addRelayCell(&relay, &Ks_cr, &Kt_cr);


    //Compute:-----------------------------------------------------------------
    S->computeSpatiotemporal();
    S->computeFourierTransform();
    io.writeStimulus(S.get());
    S->clearSpatioTemporal();

    vector<Neuron *> neurons;
//    neurons.push_back(&ganglion);
    neurons.push_back(&relay);
    neurons.push_back(&cortical);

    io.writeIntegratorProperties(&integrator);
    for(Neuron* neuron : neurons){

        neuron->computeResponse(S.get());
        io.writeResponse(neuron);
        neuron->clearResponse();

        neuron->computeImpulseResponse();
        io.writeImpulseResponse(neuron);
        neuron->clearImpulseResponse();

    }



    t = clock() - t;
    printf ("%f seconds.\n",((float)t)/CLOCKS_PER_SEC);

    return 0;
}
