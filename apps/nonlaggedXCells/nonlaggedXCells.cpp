#include <iostream>
#include <unistd.h>
#include <time.h>
#include <yaml-cpp/yaml.h>

#include <lgnSimulator.h>

using namespace std;
using namespace lgnSimulator;

int main(int argc, char* argv[])
{

    cout << "========= LGN Simulator: Nonlagged X-Cells =========" << endl;
    clock_t t;
    t = clock();

    if(argc < 2) {
        cerr << "Too few arguments." << endl;
        return 1;
    }

    //read config file-------------------------------------------------------
    YAML::Node cfg = YAML::LoadFile(argv[1]);
    const YAML::Node& ganglionImpRes = cfg["ganglionImpRes"];
    const YAML::Node& spatialKernelSettings = cfg["spatialKernels"];
    const YAML::Node& temporalKernelSettings = cfg["temporalKernels"];
    string outputFilename = cfg["outputFile"].as<std::string>();

    //Output manager:---------------------------------------------------------
    OutputManager io(&outputFilename);

    //Integrator--------------------------------------------------------------
    Integrator integrator = createIntegrator(cfg);

    //Stim---------------------------------------------------------------------
    unique_ptr<Grating> S = createGratingStimulus(&integrator, &cfg);


    //Kernels:---------------------------------------------------------
    SpatialDelta Ks_rg = createSpatialDeltaKernel(&spatialKernelSettings);
    TemporallyConstant Kt_rg = createTemporallyConstantKernel(&temporalKernelSettings);
    SeparableKernel Krg(&Ks_rg, &Kt_rg);

    DOG Ks_rig = createSpatialDOGKernel(&spatialKernelSettings);
    TemporallyConstant Kt_rig = createTemporallyConstantKernel(&temporalKernelSettings);
    SeparableKernel Krig(&Ks_rig, &Kt_rig);


    //Ganglion cell:-----------------------------------------------------------
    DOG Wg_s = createSpatialDOGKernel(&ganglionImpRes);
    TemporallyConstant Wg_t = createTemporallyConstantKernel(&ganglionImpRes);
    SeparableKernel Wg(&Wg_s, &Wg_t);
    GanglionCell ganglion(integrator, &Wg);

    //Relay cell: -------------------------------------------------------------
    RelayCell relay(integrator);

    //Connect neurons:---------------------------------------------------------
    relay.addGanglionCell(&ganglion, &Krg);
    relay.addGanglionCell(&ganglion, &Krig);

    //Compute:-----------------------------------------------------------------
    S->computeSpatiotemporal();
    S->computeFourierTransform();
    io.writeStimulus(S.get());
    S->clearSpatioTemporal();

    vector<Neuron *> neurons;
    neurons.push_back(&ganglion);
    neurons.push_back(&relay);

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

