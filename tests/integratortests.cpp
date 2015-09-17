#include <unittest++/UnitTest++.h>
#include <armadillo>
#include <iostream>


#include "stimuli/patchgrating.h"


using namespace std;
using namespace arma;

double expFunction(double x, double y, double z)
{
    return exp(-x*x - y*y - z*z);
}

double cosFunction(double x, double y, double z)
{
    return cos(2*x - y - z);
}

SUITE(integrator){

//    TEST(stimuli){
//        Config cfg;
//        cfg.readFile("../../eDOG/tests/configTests.cfg");
//        PatchGrating stim(&cfg);

//        const Setting & root = cfg.getRoot();
//        double wpg = root["stimuliSettings"]["w"];

//        int N = 1e3;
//        double *w = new double [N];
//        double *x = new double [N];
//        gauleg(-300, 300, x, w, N);


//        double I = 0;
//        double rx = 0.5;
//        double ry = 1.0;
//        double t = 3.2;
//        double Ie = stim.spatial({rx, ry}, t);

//        for(int i = 0; i < N; i++){
//            for(int j = 0; j < N; j++){
////                cout << x[i] <<  "    "<< x[j] << endl;
//                double s = stim.frequency({x[i], x[j]}, wpg);
//                I+= s * cos(rx*x[i]+ ry*x[j] - wpg * t) * w[i] * w[j];
//            }
//        }
//        I /= 8*PI*PI*PI;


//        CHECK_CLOSE(I, Ie, 1);
////        cout << setprecision(10) << Ie << '\n';
////        cout << setprecision(10) << I << '\n';


//    }




}


