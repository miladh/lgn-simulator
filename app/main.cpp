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

    //----------------------------------------------------------------------------

    PatchGrating S(&cfg);
    OutputManager io(&cfg);

    RelayCell relay(&cfg, &S);
    GanglionCell ganglion(&cfg, &S);
    CorticalCell cortical(&cfg, &S);

    DOG dog(dogA, doga, dogB, dogb);
    Gaussian gauss;
    DecayingExponential Ktg(tau_rg, 0);
    DecayingExponential Ktc(tau_rc, delay);

    relay.addGanglionCell(&ganglion,&dog, &Ktg);
    relay.addCorticalNeuron(&cortical, &gauss, &Ktc);





    double t = 0.0;
    for (int i = 0; i < nSteps; i++){
        relay.computeResponse(t);
        io.writeResponse(i, relay, S);
        cout <<"timestep: " << i << " of " << nSteps << endl;
        t+=dt;
    }




    return 0;
}
