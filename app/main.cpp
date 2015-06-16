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

    ImpulseResponse G;
    Stimuli S;
    Trapezoidal* I = new Trapezoidal((domain(0), domain(1), domain(2)));
    Response R(G, S, I, realGrid, complexGrid);

    OutputManager io(&cfg);


    for (int i = 0; i < nSteps; i++){
        io.writeResponse(i,R);
    }

    cout << R.complex() << endl;


    return 0;
}
