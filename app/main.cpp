#include <iostream>

#include <impulseResponse.h>
#include <trapezoidal.h>
#include <response.h>
#include <outputmanager.h>
#include <unistd.h>

#include "ganglion/ganglion.h"
#include "stimuli/patchgrating.h"
using namespace std;

int main()
{

    cout <<"=====Extended-DOG Model====="<<endl;

    //read config file---------------------------------------------------------------
    Config cfg;
    cfg.readFile("../../eDOG/app/config.cfg");
    const Setting & root = cfg.getRoot();

    int nSteps = root["dynamicSettings"]["nSteps"];
    double dt = root["dynamicSettings"]["dt"];

    vec mesh = zeros<vec>(3);
    vec domain = zeros<vec>(3);

    const Setting &grid = root["gridSettings"]["grid"];
    const Setting &integrationDomain = root["gridSettings"]["integrationDomain"];

    for(int i =0; i < 3; i++){
        mesh[i] = grid[i];
        domain[i] = integrationDomain[i];
    }

    //-------------------------------------------------------------------------------

    ImpulseResponse G(&cfg);
    PatchGrating S(&cfg);
    Trapezoidal* I = new Trapezoidal((domain(0), domain(1), domain(2)));
    Response R(&G, &S, I, mesh, domain);

    OutputManager io(&cfg);
    double t = 0.0;
    for (int i = 0; i < nSteps; i++){
        R.compute(t);
        io.writeResponse(i, R, G, S);
        cout <<"timestep: " << i << " of " << nSteps << endl;
        t+=dt;
    }



    return 0;
}
