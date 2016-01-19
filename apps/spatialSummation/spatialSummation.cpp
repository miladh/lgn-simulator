#include <iostream>
#include <unistd.h>
#include <time.h>
#include <yaml-cpp/yaml.h>

#include <edog.h>

using namespace std;
using namespace edog;

int main(int argc, char* argv[])
{

    cout << "=====Extended-DOG Model=====" << endl;
    clock_t t;
    t = clock();


//    string outputPath = "../../../eDOG/apps/spatialSummation/spatialSummationConfig.yaml";
    if(argc < 2) {
        cerr << "Too few arguments." << endl;
        return 1;
    }
//    if(argc >= 3) {
//        string outputPath = argv[2];
//        outputPath += "/";
//        cout << "Output will be written to " << outputPath << endl;
//    }

    //read config file-------------------------------------------------------
    YAML::Node cfg = YAML::LoadFile(argv[1]);

    //Output manager:----------------------------------------------------------
    OutputManager io(&cfg);

    //Integrator-------------------------------------------------------------
    Integrator integrator = createIntegrator(&cfg);

    //Stim---------------------------------------------------------------------
    unique_ptr<Grating> S = createGratingStimulus(&integrator, &cfg);



    //Spatial kernels:----------------------------------------------------------
    DOG dog = createDOGSpatialKernel(&cfg);


    //Temporal kernels:-------------------------------------------------------
    TemporalDelta Kt_cr(0.4);

    //Neurons:-----------------------------------------------------------------
    GanglionCell ganglion(&integrator, &dog, &Kt_cr/*, &staticNonlinearity*/);

    vector<Neuron *> neurons;
    neurons.push_back(&ganglion);

    //Compute:-----------------------------------------------------------------
    S->computeSpatiotemporal();
    S->computeFourierTransform();
    io.writeStimulus(S.get());
    S->clearSpatioTemporal();

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

