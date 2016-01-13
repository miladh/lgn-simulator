#include <iostream>

#include <outputmanager.h>
#include <unistd.h>
#include <time.h>

#include "stimuli/grating/grating.h"
#include "stimuli/staticimage.h"
#include "stimuli/naturalscenevideo.h"

#include "integrator.h"

#include "neurons/ganglioncell.h"
#include "neurons/interneuron.h"
#include "neurons/relaycell.h"
#include "neurons/corticalcell.h"

#include "spatialKernels/dog.h"
#include "spatialKernels/ellipticgaussian.h"

#include "temporalKernels/decayingexponential.h"
#include "temporalKernels/dampedoscillator.h"
#include "temporalKernels/temporallyconstant.h"
#include "temporalKernels/temporaldelta.h"

#include "staticNonlinearity/thresholdnonlinearity.h"
#include "staticNonlinearity/heavisidenonlinearity.h"
#include "staticNonlinearity/sigmoidalnonlinearity.h"



using namespace std;

int main()
{

    cout << "=====Extended-DOG Model=====" << endl;
    clock_t t;
    t = clock();


    //read config file-------------------------------------------------------
    Config cfg;
    cfg.readFile("../../eDOG/app/config.yaml");

    //Output manager:----------------------------------------------------------
    OutputManager io(&cfg);

    //Integrator-------------------------------------------------------------
    Integrator integrator = createIntegrator(&cfg);

//    Stim---------------------------------------------------------------------
//        unique_ptr<NaturalSceneVideo> S = createNaturalSceneVideoStimulus(&integrator,&cfg);
//        unique_ptr<StaticImage> S = createStaticImageStimulus(&integrator,&cfg);
        unique_ptr<Grating> S = createGratingStimulus(&integrator,&cfg);



    //Spatial kernels:----------------------------------------------------------
    DOG dog = createDOGSpatialKernel(&cfg);
//    DOG gauss(0.3, 0.05, 0, 1.);
//    EllipticGaussian ellipticGauss = createEllipticGaussianSpatialKernel(&cfg);


    //Temporal kernels:-------------------------------------------------------
//    DecayingExponential Kt_rg(0.21,0);
//    DecayingExponential Kt_rig(0.30,0);
//    DecayingExponential Kt_rc(0.42,0.10);
    TemporalDelta Kt_cr(0.4);
//    DampedOscillator damped = createDampedOscillatorTemporalKernel(&cfg);
    //    TemporallyConstant tempConst = createTemporallyConstantTemporalKernel(&cfg);


    //Static nonlinearity-------------------------------------------------------------
    ThresholdNonlinearity staticNonlinearity  = createThresholdNonlinearity(&cfg);
//    SigmoidalNonlinearity staticNonlinearity  = createSigmoidalNonlinearity(&cfg);
//    HeavisideNonlinearity staticNonlinearity;


    //Neurons:-----------------------------------------------------------------
    GanglionCell ganglion(&integrator, &dog, &Kt_cr/*, &staticNonlinearity*/);
    RelayCell relay(&integrator);
//    CorticalCell cortical(&integrator);
//    Interneuron interneuron(&integrator);

    vector<Neuron *> neurons;
    neurons.push_back(&ganglion);
    neurons.push_back(&relay);
//    neurons.push_back(&cortical);
//    neurons.push_back(&interneuron);

    //connect neurons----------------------------------------------------------
    relay.addGanglionCell(&ganglion, &dog, &Kt_cr);
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

