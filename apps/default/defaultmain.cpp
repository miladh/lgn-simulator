#include <iostream>
#include <yaml-cpp/yaml.h>
#include <unistd.h>
#include <time.h>

#include <lgnSimulator.h>

using namespace std;
using namespace lgnSimulator;

int main()
{

    cout << "=====LGN Simulator Model=====" << endl;
    clock_t t;
    t = clock();


    //read config file-------------------------------------------------------
    YAML::Node cfg = YAML::LoadFile("../../../lgn-simulator/apps/default/defaultConfig.yaml");

    //Output manager:----------------------------------------------------------
    string outputFilename = cfg["outputFile"].as<std::string>();
    OutputManager io(&outputFilename);

    //Integrator-------------------------------------------------------------
    Integrator integrator = createIntegrator(cfg);

//    Stim---------------------------------------------------------------------
//        unique_ptr<NaturalSceneVideo> S = createNaturalSceneVideoStimulus(&integrator,&cfg);
//        unique_ptr<StaticImage> S = createStaticImageStimulus(&integrator,&cfg);
        unique_ptr<Grating> S = createGratingStimulus(integrator, cfg);



    //Spatial kernels:----------------------------------------------------------
    DOG dog = createSpatialDOGKernel(cfg);
//    DOG gauss(0.3, 0.05, 0, 1.);
//    EllipticGaussian ellipticGauss = createEllipticGaussianSpatialKernel(&cfg);


    //Temporal kernels:-------------------------------------------------------
//    DecayingExponential Kt_rg(0.21,0);
//    DecayingExponential Kt_rig(0.30,0);
//    DecayingExponential Kt_rc(0.42,0.10);
//    DampedOscillator damped = createDampedOscillatorTemporalKernel(&cfg);
//        TemporallyConstant tempConst = createTemporallyConstantTemporalKernel(&cfg);
    TemporalDelta Kt_cr(0.4);

    //Static nonlinearity-------------------------------------------------------------
//    ThresholdNonlinearity staticNonlinearity  = createThresholdNonlinearity(&cfg);
//    SigmoidalNonlinearity staticNonlinearity  = createSigmoidalNonlinearity(&cfg);
//    HeavisideNonlinearity staticNonlinearity;


    //Kernels:-------------------------------------------------------
    SeparableKernel Kt(&dog, &Kt_cr);

    //Neurons:-----------------------------------------------------------------
    GanglionCell ganglion(integrator, &Kt/*, &staticNonlinearity*/);
    RelayCell relay(integrator);
//    CorticalCell cortical(&integrator);
//    Interneuron interneuron(&integrator);

    vector<Neuron *> neurons;
    neurons.push_back(&ganglion);
    neurons.push_back(&relay);
//    neurons.push_back(&cortical);
//    neurons.push_back(&interneuron);

    //connect neurons----------------------------------------------------------
    relay.addGanglionCell(&ganglion, &Kt);
//    relay.addCorticalNeuron(&cortical, &ellipticGauss, &Kt_rc);
//    relay.addInterNeuron(&interneuron,&gauss, &damped);

//    interneuron.addGanglionCell(&ganglion, &gauss, &Kt_rig);
//    interneuron.addCorticalNeuron(&cortical, &ellipticGauss, &Kt_rc);

//    cortical.addRelayCell(&relay, &gauss, &Kt_cr);



    //Compute:-----------------------------------------------------------------
    S->computeSpatiotemporal();
    S->computeFourierTransform();
    io.writeStimulus(S.get());
    S->clearSpatioTemporal();

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

