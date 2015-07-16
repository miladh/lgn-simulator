#include <iostream>

#include <trapezoidal.h>
#include <outputmanager.h>
#include <unistd.h>

#include "stimuli/patchgrating.h"

#include "neurons/relaycell.h"
#include "neurons/ganglioncell.h"
#include "neurons/corticalcell.h"

#include "spatialKernels/dog.h"
#include "spatialKernels/gaussian.h"

#include "temporalKernels/decayingexponential.h"
#include "temporalKernels/diracDelta.h"


using namespace std;

int main()
{

    cout << "=====Extended-DOG Model=====" << endl;

    //read config file---------------------------------------------------------------
    Config cfg;
    cfg.readFile("../../eDOG/app/config.cfg");
    const Setting & root = cfg.getRoot();

    int nSteps = root["dynamicSettings"]["nSteps"];
    double dt = root["dynamicSettings"]["dt"];

    double dogA = root["dogSettings"]["A"];
    double doga = root["dogSettings"]["a"];
    double dogB = root["dogSettings"]["B"];
    double dogb = root["dogSettings"]["b"];
    double tau_rg = root["temporalSettings"]["tau_rg"];
    double tau_rc = root["temporalSettings"]["tau_rc"];
    double delay = root["temporalSettings"]["delay"];
    double weight = root["loopKernelSettings"]["C"];
    double spread = root["loopKernelSettings"]["c"];


    //----------------------------------------------------------------------------

    PatchGrating S(&cfg);
    OutputManager io(&cfg);

    //Kernels:
    DOG dog(dogA, doga, dogB, dogb);
    Gaussian gauss(weight, spread);
    DecayingExponential Ktg(tau_rg, 0);
    DecayingExponential Ktc(tau_rc, delay);
    DiracDelta delta(0.0);


    //Neurons:
    RelayCell relay(&cfg, &S);
    CorticalCell cortical(&cfg, &S);
    GanglionCell ganglion(&cfg, &S, &dog, &Ktg);

    vector<Neuron *> neurons;
    neurons.push_back(&relay);
    neurons.push_back(&cortical);



    relay.addGanglionCell(&ganglion,&dog, &delta);
    relay.addCorticalNeuron(&cortical, &gauss, &Ktc);

    cortical.addRelayCell(&relay,&dog, &delta);


    double t = 0.0;
    for (int i = 0; i < nSteps; i++){
        cortical.computeResponse(t);
        relay.computeResponse(t);
        io.writeResponse(i, neurons, S);
        cout <<"timestep: " << i << " of " << nSteps << endl;
        t+=dt;
    }




    return 0;
}
