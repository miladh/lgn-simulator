#include <iostream>

#include <outputmanager.h>
#include <unistd.h>

#include "stimuli/patchgrating.h"
#include "stimuli/grating.h"
#include "stimuli/oscillatinggaussian.h"

#include "integrator.h"

#include "neurons/ganglioncell.h"
#include "neurons/interneuron.h"
#include "neurons/relaycell.h"
#include "neurons/corticalcell.h"

#include "spatialKernels/dog.h"
#include "spatialKernels/gaussian.h"
#include "spatialKernels/ellipticgaussian.h"

#include "temporalKernels/decayingexponential.h"
#include "temporalKernels/dampedoscillator.h"


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
//    Grating S = createGratingStimulus(&integrator,&cfg);
    PatchGrating S = createPatchGratingStimulus(&integrator,&cfg);
//    OscillatingGaussian S = createOscillatingGaussianStimulus(&integrator,&cfg);



    //Spatial kernels:----------------------------------------------------------
    DOG dog = createDOGSpatialKernel(&cfg);
    Gaussian gauss = createGaussianSpatialKernel(&cfg);
    EllipticGaussian ellipticGauss = createEllipticGaussianSpatialKernel(&cfg);


    //Temporal kernels:-------------------------------------------------------
    DecayingExponential Kt = createDecayingExponentialTemporalKernel(&cfg);
    DampedOscillator damped = createDampedOscillatorTemporalKernel(&cfg);


    //Neurons:-----------------------------------------------------------------
    GanglionCell ganglion(&integrator, &dog, &damped);
    RelayCell relay(&integrator);
    Interneuron interneuron(&integrator);
    CorticalCell cortical(&integrator);

    vector<Neuron *> neurons;
    neurons.push_back(&ganglion);
    neurons.push_back(&relay);
    neurons.push_back(&interneuron);
    neurons.push_back(&cortical);



    //connect neurons----------------------------------------------------------
    relay.addGanglionCell(&ganglion,&gauss, &damped);
    relay.addInterNeuron(&interneuron,&dog, &damped);
    relay.addCorticalNeuron(&cortical, &ellipticGauss, &Kt);

    interneuron.addGanglionCell(&ganglion,&dog, &Kt);
    interneuron.addCorticalNeuron(&cortical, &ellipticGauss, &Kt);

    cortical.addRelayCell(&relay, &dog, &Kt);



    //Compute:-----------------------------------------------------------------
    S.computeSpatiotemporal();
    S.computeFourierTransform();

    ganglion.computeResponse(&S);
    ganglion.computeImpulseResponse();

    relay.computeResponse(&S);
    relay.computeImpulseResponse();

    interneuron.computeResponse(&S);
    interneuron.computeImpulseResponse();

    cortical.computeResponse(&S);
    cortical.computeImpulseResponse();





    //Output manager:----------------------------------------------------------
    OutputManager io(&cfg);
    io.writeResponse(neurons, S);


    return 0;
}

