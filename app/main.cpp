#include <iostream>

#include <impulseResponse.h>
#include <stimuli.h>
#include <trapezoidal.h>
#include <response.h>
#include <outputmanager.h>
#include <unistd.h>
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

    vec realGrid = zeros<vec>(3);
    vec complexGrid = zeros<vec>(3);
    vec domain = zeros<vec>(3);

    const Setting &real = root["gridSettings"]["realGrid"];
    const Setting &complex = root["gridSettings"]["complexGrid"];
    const Setting &integrationDomain = root["gridSettings"]["integrationDomain"];

    for(int i =0; i < 3; i++){
        realGrid[i] = real[i];
        complexGrid[i] = complex[i];
        domain[i] = integrationDomain[i];
    }

    //-------------------------------------------------------------------------------

    ImpulseResponse G(&cfg);
    Stimuli S(&cfg);
    Trapezoidal* I = new Trapezoidal((domain(0), domain(1), domain(2)));
    Response R(G, S, I, realGrid, complexGrid);

    OutputManager io(&cfg);

    double t = 0.0;
    for (int i = 0; i < nSteps; i++){
        R.computeComplex(t);
        io.writeResponse(i,R);
        cout << R.complex() << endl;
        t+=dt;
    }



    return 0;
}
