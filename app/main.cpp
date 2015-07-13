#include <iostream>

#include <trapezoidal.h>
#include <outputmanager.h>
#include <unistd.h>

#include "ganglion/ganglion.h"
#include "stimuli/patchgrating.h"
#include "relay/originaledog.h"

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


    //----------------------------------------------------------------------------


    DOG dog(1., 1., 0.1, 0.2);
    GanglionDOG ganglion(&dog);

    PatchGrating S(&cfg);
    OriginalEDOG R(&cfg, &ganglion, &S);

    OutputManager io(&cfg);
    double t = 0.0;
    for (int i = 0; i < nSteps; i++){
        R.computeResponse(t);
//        io.writeResponse(i, R, G, S);
        cout <<"timestep: " << i << " of " << nSteps << endl;
        t+=dt;
    }




    return 0;
}
