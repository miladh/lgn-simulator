#include <iostream>

#include <outputmanager.h>
#include <unistd.h>

#include "stimuli/patchgrating.h"<l
#include "stimuli/grating.h"

#include "integrator.h"

#include "neurons/ganglioncell.h"
#include "neurons/interneuron.h"
#include "neurons/relaycell.h"
#include "neurons/corticalcell.h"

#include "spatialKernels/dog.h"
#include "spatialKernels/gaussian.h"
#include "spatialKernels/ellipticgaussian.h"

#include "temporalKernels/decayingexponential.h"
#include "temporalKernels/diracDelta.h"
#include "temporalKernels/dampedoscillator.h"


using namespace std;

int main()
{

    cout << "=====Extended-DOG Model=====" << endl;

    //read config file---------------------------------------------------------------
    Config cfg;
    cfg.readFile("../../eDOG/app/config.cfg");
    const Setting & root = cfg.getRoot();


    double dogA = root["dogSettings"]["A"];
    double doga = root["dogSettings"]["a"];
    double dogB = root["dogSettings"]["B"];
    double dogb = root["dogSettings"]["b"];
    double tau_rg = root["temporalSettings"]["tau_rg"];
    double tau_rc = root["temporalSettings"]["tau_rc"];
    double delay = root["temporalSettings"]["delay"];
    double weight = 0.5;
    double spread = 1.1;

    double contrast = root["stimuliSettings"]["C"];
    double spotDiameter = root["stimuliSettings"]["d"];


    //----------------------------------------------------------------------------
    Integrator integrator = createIntegrator(&cfg);
    OutputManager io(&cfg);


    //Stim---------------------------------------------------------------------
//    Grating S = createGratingStimulus(&integrator,&cfg);
    PatchGrating S = createPatchGratingStimulus(&integrator,&cfg);


    //Spatial kernels:----------------------------------------------------------
    DOG dog(dogA, doga, dogB, dogb);
    Gaussian gauss(weight, spread);
    EllipticGaussian ellipticGauss(weight, PI/4*3, 0.1, 1.0);


    //Temporal kernels:-------------------------------------------------------
    DecayingExponential Ktg(tau_rg, 0);
    DecayingExponential Ktc(tau_rc, delay);
    DiracDelta delta(0.0);
    DampedOscillator damped(10, 1.38);

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
    relay.addCorticalNeuron(&cortical, &ellipticGauss, &Ktc);

    interneuron.addGanglionCell(&ganglion,&dog, &Ktg);
    interneuron.addCorticalNeuron(&cortical, &ellipticGauss, &Ktc);

    cortical.addRelayCell(&relay, &dog, &Ktg);






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





    //write:-------------------------------------------------------------------
    io.writeResponse(neurons, S);


    return 0;
}

