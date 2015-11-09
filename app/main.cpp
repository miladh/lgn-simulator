#include <iostream>

#include <outputmanager.h>
#include <unistd.h>

#include "stimuli/patchgrating.h"
#include "stimuli/grating.h"
#include "stimuli/oscillatinggaussian.h"
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


using namespace std;

int main()
{

    cout << "=====Extended-DOG Model=====" << endl;

    //read config file-------------------------------------------------------
    Config cfg;
    cfg.readFile("../../eDOG/app/config.cfg");

    //Integrator-------------------------------------------------------------
    Integrator integrator = createIntegrator(&cfg);

    //Stim---------------------------------------------------------------------
//    NaturalSceneVideo S = createNaturalSceneVideoStimulus(&integrator,&cfg);
    StaticImage S = createStaticImageStimulus(&integrator,&cfg);
//    Grating S = createGratingStimulus(&integrator,&cfg);
//    PatchGrating S = createPatchGratingStimulus(&integrator,&cfg);
//    OscillatingGaussian S = createOscillatingGaussianStimulus(&integrator,&cfg);



    //Spatial kernels:----------------------------------------------------------
    DOG dog = createDOGSpatialKernel(&cfg);
    EllipticGaussian ellipticGauss = createEllipticGaussianSpatialKernel(&cfg);


    //Temporal kernels:-------------------------------------------------------
    DecayingExponential Kt = createDecayingExponentialTemporalKernel(&cfg);
    DampedOscillator damped = createDampedOscillatorTemporalKernel(&cfg);
    TemporallyConstant tempConst = createTemporallyConstantTemporalKernel(&cfg);


    //Neurons:-----------------------------------------------------------------
    GanglionCell ganglion(&integrator, &dog, &Kt);
    RelayCell relay(&integrator);
//    Interneuron interneuron(&integrator);
//    CorticalCell cortical(&integrator);

    vector<Neuron *> neurons;
    neurons.push_back(&ganglion);
    neurons.push_back(&relay);
//    neurons.push_back(&interneuron);
//    neurons.push_back(&cortical);



    //connect neurons----------------------------------------------------------
    relay.addGanglionCell(&ganglion,&dog, &Kt);
//    relay.addInterNeuron(&interneuron,&gauss, &damped);
//    relay.addCorticalNeuron(&cortical, &ellipticGauss, &Kt);

//    interneuron.addGanglionCell(&ganglion,&gauss, &Kt);
//    interneuron.addCorticalNeuron(&cortical, &ellipticGauss, &Kt);

//    cortical.addRelayCell(&relay, &dog, &damped);



    //Compute:-----------------------------------------------------------------
    S.computeSpatiotemporal();
    S.computeFourierTransform();

    ganglion.computeResponse(&S);
    relay.computeResponse(&S);
//    interneuron.computeResponse(&S);
//    cortical.computeResponse(&S);





    //Output manager:----------------------------------------------------------
    OutputManager io(&cfg);
    io.writeResponse(neurons, S);


    return 0;
}

