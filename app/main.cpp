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
    vec realGrid = {-4, 4, 10};
    vec fourierGrid = {-3, 3, 10};
    vec domain = {-1, 1, 10};


    ImpulseResponse G;
    Stimuli S;
    Trapezoidal* I = new Trapezoidal((domain(0), domain(1), domain(2)));
    Response R(G, S, I, realGrid, fourierGrid);


    Config cfg;
    OutputManager io(&cfg);


    for (int i = 0; i < 2; i++){
        io.writeResponse(i,R);
    }



    cout << R.complex() << endl;


    return 0;
}
