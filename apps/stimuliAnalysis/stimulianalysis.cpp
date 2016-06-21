#include <iostream>
#include <unistd.h>
#include <time.h>
#include <yaml-cpp/yaml.h>

#include <lgnSimulator.h>

using namespace std;
using namespace lgnSimulator;

int main(int argc, char* argv[])
{

    cout << "=====LGN Simulator Model: Stimuli Analysis =====" << endl;
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
    io.copyConfigFile(argv[1]);

    //Integrator--------------------------------------------------------------
    Integrator integrator = createIntegrator(cfg["grid"]);
    io.writeIntegratorProperties(integrator);

    //Stim---------------------------------------------------------------------
    unique_ptr<Grating> S = createGratingStimulus(integrator, cfg["stimulus"]);


    //Compute:-----------------------------------------------------------------
    io.writeStimulusProperties(S.get());

    S->computeFourierTransform();
    S->computeSpatiotemporal();
    io.writeStimulus(S.get());


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

