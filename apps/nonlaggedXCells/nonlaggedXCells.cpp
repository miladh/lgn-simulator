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


    //Output manager:---------------------------------------------------------
    OutputManager io(&cfg);

    //Integrator--------------------------------------------------------------
    Integrator integrator = createIntegrator(&cfg);

    //Stim---------------------------------------------------------------------
    unique_ptr<Grating> S = createGratingStimulus(&integrator, &cfg);


    //Spatial kernels:---------------------------------------------------------
    SpatialDelta Ks_rg = createSpatialDeltaSpatialKernel(&cfg);
    SpatiallyConstant Ks_ri = createSpatiallyConstantSpatialKernel(&cfg);
    DOG Ks_ig = createDOGSpatialKernel(&cfg);

    //Temporal kernels:--------------------------------------------------------
    TemporalDelta Kt_rg = createTemporalDeltaKernel(&cfg);
    TemporalDelta Kt_ri = createTemporalDeltaKernel(&cfg);
    TemporalDelta Kt_ig = createTemporalDeltaKernel(&cfg);


    //Ganglion cell:-----------------------------------------------------------
    DOG Wg_s = createDOGSpatialKernel(&cfg);
    TemporalDelta Wg_t = createTemporalDeltaKernel(&cfg);
    GanglionCell ganglion(&integrator, &Wg_s, &Wg_t);

    //Relay cell: -------------------------------------------------------------
    RelayCell relay(&integrator);

    //Interneuron:-------------------------------------------------------------
    Interneuron interneuron(&integrator);


    //Connect neurons:---------------------------------------------------------
    relay.addGanglionCell(&ganglion, &Ks_rg, &Kt_rg);
    relay.addInterNeuron(&interneuron, &Ks_ri, &Kt_ri);
    interneuron.addGanglionCell(&ganglion, &Ks_ig, &Kt_ig);

    //Compute:-----------------------------------------------------------------
    S->computeSpatiotemporal();
    S->computeFourierTransform();
    io.writeStimulus(S.get());
    S->clearSpatioTemporal();

    vector<Neuron *> neurons;
    neurons.push_back(&ganglion);
    neurons.push_back(&relay);
    neurons.push_back(&interneuron);

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

