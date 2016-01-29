#include <iostream>
#include <unistd.h>
#include <time.h>
#include <yaml-cpp/yaml.h>

#include <edog.h>

using namespace std;
using namespace edog;

int main(int argc, char* argv[])
{

    cout << "========= Extended-DOG: Phase-reversed push-pull =========" << endl;
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

    //Output manager:---------------------------------------------------------
    OutputManager io(&cfg);

    //Integrator--------------------------------------------------------------
    Integrator integrator = createIntegrator(&cfg);

    //Stim---------------------------------------------------------------------
    unique_ptr<Grating> S = createGratingStimulus(&integrator, &cfg);

    //Ganglion cell:-----------------------------------------------------------
    DOG Wg_s = createDOGSpatialKernel(&ganglionImpRes);
    TemporallyConstant Wg_t = createTemporallyConstantKernel(&ganglionImpRes);
    GanglionCell ganglion(&integrator, &Wg_s, &Wg_t);

    //Relay cell: -------------------------------------------------------------
    RelayCell relay(&integrator);
    CorticalCell cortical(&integrator);

    //Spatial kernels:---------------------------------------------------------
    SpatialDelta Ks_rg = createSpatialDeltaKernel(&Ks_rgSettings);
    SpatialDelta Ks_cr = createSpatialDeltaKernel(&Ks_crSettings);
    DOG Ks_rc = createDOGSpatialKernel(&Ks_rcSettings);

    //Temporal kernels:--------------------------------------------------------
    TemporallyConstant Kt_rg = createTemporallyConstantKernel(&Kt_rgSettings);
    TemporallyConstant Kt_cr = createTemporallyConstantKernel(&Kt_crSettings);
    TemporalDelta Kt_rc = createTemporalDeltaKernel(&Kt_rcSettings);


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
    neurons.push_back(&ganglion);
    neurons.push_back(&relay);
    neurons.push_back(&cortical);

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

