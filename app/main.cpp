#include <iostream>

#include <impulseResponse.h>
#include <stimuli.h>
#include <trapezoidal.h>
#include <response.h>

using namespace std;

int main()
{


    cout <<"=====Extended-DOG Model====="<<endl;
    ImpulseResponse G;
    Stimuli S;
    Trapezoidal* I;

    vec grid = {-1, 1, 10};

    Response R(G, S, I, grid);

    cout << R.response() << endl;




    return 0;
}
