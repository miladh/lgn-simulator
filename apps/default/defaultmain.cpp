#include <iostream>
#include <unistd.h>
#include <time.h>
#include <yaml-cpp/yaml.h>

#include <lgnSimulator.h>

using namespace std;
using namespace lgnSimulator;

int main(int argc, char* argv[]){

    cout << "========= LGN Simulator: Default =========" << endl;
    clock_t t;
    t = clock();

    if(argc < 2) {
        cerr << "Too few arguments." << endl;
        return 1;
    }

    //read config file-------------------------------------------------------
    YAML::Node cfg = YAML::LoadFile(argv[1]);
    const YAML::Node& ganglionImpRes = cfg["ganglionImpRes"];
    const YAML::Node& Ks_rgSettings = cfg["kernels"]["Krg"]["spatial"];
    const YAML::Node& Kt_rgSettings = cfg["kernels"]["Krg"]["temporal"];
    const YAML::Node& Ks_rcSettings = cfg["kernels"]["Krc"]["spatial"];
    const YAML::Node& Kt_rcSettings = cfg["kernels"]["Krc"]["temporal"];
    const YAML::Node& K_crSettings = cfg["kernels"]["Kcr"];


    string outputFilename = cfg["outputFile"].as<std::string>();
    double Rg_0 = cfg["backgroundResponse"]["Rg_0"].as<double>();
    double Rr_0 = cfg["backgroundResponse"]["Rr_0"].as<double>();
    double Rc_0 = cfg["backgroundResponse"]["Rc_0"].as<double>();

    //Output manager:---------------------------------------------------------
    OutputManager io(outputFilename);

    //Integrator--------------------------------------------------------------
    Integrator integrator = createIntegrator(cfg);

    //Stim---------------------------------------------------------------------
    unique_ptr<Grating> S = createGratingStimulus(&integrator, cfg);

//    //Ganglion cell:-----------------------------------------------------------
//    DOG Wg_s = createSpatialDOGKernel(ganglionImpRes);
//    Biphasic Wg_t = createTemporalBiphasicKernel(ganglionImpRes);

//    SeparableKernel Wg(&Wg_s, &Wg_t);
//    GanglionCell ganglion(&integrator, Wg, Rg_0);

//    //Relay cell: -------------------------------------------------------------
//    RelayCell relay(&integrator, Rr_0);
//    CorticalCell cortical(&integrator, Rc_0);

//    //Kernels:---------------------------------------------------------
//    SpatialDelta Ks_rg = createSpatialDeltaKernel(Ks_rgSettings);
//    TemporalDelta Kt_rg = createTemporalDeltaKernel(Kt_rgSettings);
//    SeparableKernel Krg(&Ks_rg, &Kt_rg);

//    SpatialGaussian Ks_rc = createSpatialGaussianKernel(Ks_rcSettings);
//    DOE Kt_rc = createTemporalDOEKernel(Kt_rcSettings);
//    SeparableKernel Krc(&Ks_rc, &Kt_rc);

//    NonseparableDOG Kcr = createNonseparableDOGKernel(K_crSettings);



//    //Connect neurons:---------------------------------------------------------
//    relay.addGanglionCell(&ganglion, Krg);
//    relay.addCorticalNeuron(&cortical, Krc);
//    cortical.addRelayCell(&relay, Kcr);

//    //Compute:-----------------------------------------------------------------
//    S->computeSpatiotemporal();
//    S->computeFourierTransform();
//    io.writeStimulus(S.get());
//    S->clearSpatioTemporal();

//    vector<Neuron *> neurons;
//    neurons.push_back(&ganglion);
//    neurons.push_back(&relay);
//    neurons.push_back(&cortical);

//    io.writeIntegratorProperties(integrator);
//    for(Neuron* neuron : neurons){

//        neuron->computeResponse(S.get());
//        io.writeResponse(*neuron);
//        neuron->clearResponse();

//        neuron->computeImpulseResponse();
//        io.writeImpulseResponse(*neuron);
//        neuron->clearImpulseResponse();

//    }

//    t = clock() - t;
//    double elapsedTime = ((float)t)/CLOCKS_PER_SEC;
//    if(elapsedTime <= 60){
//        printf ("%f seconds.\n", elapsedTime);
//    }else{
//        printf ("%f minutes.\n", elapsedTime/60);
//    }

    return 0;
}

