#include <unittest++/UnitTest++.h>
#include <armadillo>
#include <iostream>


#include "stimuli/patchgrating.h"
#include <lib.h>


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

    TEST(stimuli){
        Config cfg;
        cfg.readFile("../../eDOG/tests/configTests.cfg");
        PatchGrating stim(&cfg);

        const Setting & root = cfg.getRoot();
        double wpg = root["stimuliSettings"]["w"];

        int N = 1e3;
        double *w = new double [N];
        double *x = new double [N];
        gauleg(-300, 300, x, w, N);


        double I = 0;
        double rx = 0.5;
        double ry = 1.0;
        double t = 3.2;
        double Ie = stim.real({rx, ry}, t);

        for(int i = 0; i < N; i++){
            for(int j = 0; j < N; j++){
//                cout << x[i] <<  "    "<< x[j] << endl;
                double s = stim.complex({x[i], x[j]}, wpg);
                I+= s * cos(rx*x[i]+ ry*x[j] - wpg * t) * w[i] * w[j];
            }
        }
        I /= 8*PI*PI*PI;


        CHECK_CLOSE(I, Ie, 1);
        cout << setprecision(10) << Ie << '\n';
        cout << setprecision(10) << I << '\n';


    }




    TEST(gaussLegendre) {

        int N = 100;
        double *w = new double [N];
        double *x = new double [N];
        double Iexp = 0;
        double Icos = 0;
        gauleg(-10, 10, x, w, N);


        for(int i = 0; i < N; i++){
            for(int j = 0; j < N; j++){
                for(int k = 0; k < N; k++){
                    Iexp += w[i] * w[j] * w[k] * expFunction(x[i], x[j], x[k]);
                    Icos += w[i] * w[j] * w[k] * cosFunction(x[i], x[j], x[k]);
                }
            }
        }
        CHECK_CLOSE(Iexp, 5.568327996831, 1e-11);
        CHECK_CLOSE(Icos, 1.0807773225605, 1e-11);

    }

}


