#include <iostream>

#include <outputmanager.h>
#include <unistd.h>
#include <time.h>

#include "stimuli/patchgrating.h"
#include "stimuli/grating.h"
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



using namespace std;

int main()
{

    cout << "=====Extended-DOG Model=====" << endl;
    clock_t t;
    t = clock();


    //read config file-------------------------------------------------------
    Config cfg;
    cfg.readFile("../../eDOG/app/config.cfg");

    //Output manager:----------------------------------------------------------
    OutputManager io(&cfg);

    //Integrator-------------------------------------------------------------
    Integrator integrator = createIntegrator(&cfg);

    //Stim---------------------------------------------------------------------
//        NaturalSceneVideo S = createNaturalSceneVideoStimulus(&integrator,&cfg);
    //    StaticImage S = createStaticImageStimulus(&integrator,&cfg);
        Grating S = createGratingStimulus(&integrator,&cfg);
//    PatchGrating S = createPatchGratingStimulus(&integrator,&cfg);
//        OscillatingGaussian S = createOscillatingGaussianStimulus(&integrator,&cfg);



    //Spatial kernels:----------------------------------------------------------
    DOG dog = createDOGSpatialKernel(&cfg);
    DOG gauss(1.0, 0.2, 0, 0.1);
    EllipticGaussian ellipticGauss = createEllipticGaussianSpatialKernel(&cfg);


    //Temporal kernels:-------------------------------------------------------
    DecayingExponential Kt_rg(0.21,0);
    DecayingExponential Kt_rig(0.30,0);
    DecayingExponential Kt_rc(0.42,0.10);
    TemporalDelta Kt_cr(0.0);
    DampedOscillator damped = createDampedOscillatorTemporalKernel(&cfg);
    //    TemporallyConstant tempConst = createTemporallyConstantTemporalKernel(&cfg);


    //Neurons:-----------------------------------------------------------------
    GanglionCell ganglion(&integrator, &gauss, &Kt_cr);
    RelayCell relay(&integrator);
    CorticalCell cortical(&integrator);
    Interneuron interneuron(&integrator);

    vector<Neuron *> neurons;
    neurons.push_back(&ganglion);
    neurons.push_back(&relay);
    neurons.push_back(&cortical);
    neurons.push_back(&interneuron);



    //connect neurons----------------------------------------------------------
    relay.addGanglionCell(&ganglion,&dog, &Kt_rg);
    relay.addCorticalNeuron(&cortical, &gauss, &Kt_rc);
    relay.addInterNeuron(&interneuron,&gauss, &damped);

    interneuron.addGanglionCell(&ganglion, &gauss, &Kt_rg);
    interneuron.addCorticalNeuron(&cortical, &ellipticGauss, &Kt_rc);

    cortical.addRelayCell(&relay, &dog, &Kt_cr);



    //Compute:-----------------------------------------------------------------
    S.computeSpatiotemporal();
    S.computeFourierTransform();
    io.writeStimulus(&S);
    S.clearSpatioTemporal();

    for(Neuron* neuron : neurons){

        neuron->computeResponse(&S);
        io.writeResponse(neuron);
        neuron->clearResponse();

        if(neuron->cellType() == "ganglion"){
            ganglion.GanglionCell::computeImpulseResponse();
        }else{
            neuron->computeImpulseResponse();
        }
        io.writeImpulseResponse(neuron);
        neuron->clearImpulseResponse();

    }


    t = clock() - t;
    printf ("%f seconds.\n",((float)t)/CLOCKS_PER_SEC);

    return 0;
}

