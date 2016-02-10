/**********************************************************************
 *  Test: 3d inverse fourier transform of temporal exponential function
 *
 *  Analytic source: by hand
 *
 * ********************************************************************/

#include <unittest++/UnitTest++.h>
#include <armadillo>
#include <iostream>

#include "integrator.h"

using namespace std;
using namespace arma;
using namespace lgnSimulator;

double exponential(double a, double t){

    return exp(-a* abs(t));
}

complex<double> exponentialFT(double a, double w){

    return  2. *a / (a*a + w*w);
}


SUITE(INTEGRATOR){

    TEST(exponential){
        //Mesh
        int ns = 2;
        int nt = 15;
        
        double dt = 0.001;
        double a = 0.6;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double ds = 0.01;

        Integrator integrator(nt, dt, ns, ds);

        vec s = integrator.spatialVec();
        vec k = integrator.spatialFreqVec();
        vec t = integrator.timeVec();
        vec w = integrator.temporalFreqVec();



        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = exponential(a, t[l]);
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = exponentialFT(a, w[l])
                            * Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0);
                }
            }
        }

        f *= 4.*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);
//

        // Test
//        for(int l = 0; l < Nt; l++){
//            for(int i = 0; i < Ns; i++){
//                for(int j = 0; j < Ns; j++){
//                    CHECK_CLOSE(real(g(i,j,l)),
//                                abs(G(i,j,l)), 1e-3);

//                }
//            }
//        }
    }


}
