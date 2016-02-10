#include <iostream>
#include <unistd.h>
#include <time.h>
#include <yaml-cpp/yaml.h>

#include <lgnSimulator.h>

using namespace std;
using namespace lgnSimulator;

int main(int argc, char* argv[])
{

    cout << "=====Extended-DOG Model=====" << endl;
    clock_t t;
    t = clock();

    if(argc < 2) {
        cerr << "Too few arguments." << endl;
        return 1;
    }

    //read config file-------------------------------------------------------
    YAML::Node cfg = YAML::LoadFile(argv[1]);
    string outputFilename = cfg["outputFile"].as<std::string>();


    //Output manager:----------------------------------------------------------
    OutputManager io(&outputFilename);

    //Integrator-------------------------------------------------------------
    Integrator integrator = createIntegrator(&cfg);

    //Stim---------------------------------------------------------------------
    unique_ptr<Grating> S = createGratingStimulus(&integrator, &cfg);


    //Spatial kernels:----------------------------------------------------------
    DOG dog = createDOGSpatialKernel(&cfg);


    //Temporal kernels:-------------------------------------------------------
    TemporalDelta temporalDelta = createTemporalDeltaKernel(&cfg);

    //Neurons:-----------------------------------------------------------------
    GanglionCell ganglion(&integrator, &dog, &temporalDelta);



    //Compute:-----------------------------------------------------------------
    S->computeSpatiotemporal();
    S->computeFourierTransform();
    io.writeStimulus(S.get());
    S->clearSpatioTemporal();

    ganglion.computeResponse(S.get());
    io.writeResponse(&ganglion);
    ganglion.clearResponse();

    ganglion.computeImpulseResponse();
    io.writeImpulseResponse(&ganglion);
    ganglion.clearImpulseResponse();

    io.writeIntegratorProperties(&integrator);

    t = clock() - t;
    printf ("%f seconds.\n",((float)t)/CLOCKS_PER_SEC);

    return 0;
}

