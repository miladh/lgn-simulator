#include <iostream>

#include <impulseResponse.h>
#include <stimuli.h>
#include <trapezoidal.h>
#include <response.h>

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

    cout << R.fourierSpace() << endl;


    return 0;
}
