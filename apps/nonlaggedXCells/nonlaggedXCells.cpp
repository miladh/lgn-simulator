#include <iostream>
#include <unistd.h>
#include <time.h>
#include <yaml-cpp/yaml.h>

#include <edog.h>

using namespace std;
using namespace edog;

int main(int argc, char* argv[])
{

    cout << "========= Extended-DOG: Nonlagged X-Cells =========" << endl;
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

    //Output manager:---------------------------------------------------------
    OutputManager io(&cfg);

    //Integrator--------------------------------------------------------------
    Integrator integrator = createIntegrator(&cfg);

    //Stim---------------------------------------------------------------------
    unique_ptr<Grating> S = createGratingStimulus(&integrator, &cfg);


    //Spatial kernels:---------------------------------------------------------
    SpatialDelta Ks_rg = createSpatialDeltaKernel(&spatialKernelSettings);
    DOG Ks_rig = createDOGSpatialKernel(&spatialKernelSettings);

    //Temporal kernels:--------------------------------------------------------
    TemporallyConstant Kt_rg = createTemporallyConstantKernel(&temporalKernelSettings);
    TemporallyConstant Kt_rig = createTemporallyConstantKernel(&temporalKernelSettings);


    //Ganglion cell:-----------------------------------------------------------
    DOG Wg_s = createDOGSpatialKernel(&ganglionImpRes);
    TemporallyConstant Wg_t = createTemporallyConstantKernel(&ganglionImpRes);
    GanglionCell ganglion(&integrator, &Wg_s, &Wg_t);

    //Relay cell: -------------------------------------------------------------
    RelayCell relay(&integrator);

    //Connect neurons:---------------------------------------------------------
    relay.addGanglionCell(&ganglion, &Ks_rg, &Kt_rg);
    relay.addGanglionCell(&ganglion, &Ks_rig, &Kt_rig);

    //Compute:-----------------------------------------------------------------
    S->computeSpatiotemporal();
    S->computeFourierTransform();
    io.writeStimulus(S.get());
    S->clearSpatioTemporal();

    vector<Neuron *> neurons;
    neurons.push_back(&ganglion);
    neurons.push_back(&relay);

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

