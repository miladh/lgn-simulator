#include <unittest++/UnitTest++.h>
#include <armadillo>
#include <iostream>

#include "integrator.h"
#include "math/functions.h"

#include "stimuli/oscillatinggaussian.h"

using namespace std;
using namespace arma;


double gauss(double a, vec r){
    return exp(-a*dot(r,r));
}

double gaussFT(double a, vec k){
    return PI/a*exp(-dot(k,k)/4./a);
}

SUITE(INTEGRATOR){


    TEST(gaussSpatialCosineTemporal_1){
        //Mesh
        int ns = 7;
        int nt = 4;
        double ds = 0.1;
        double dt = 0.1;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        Integrator integrator(nt, dt, ns, ds);

        vec s = integrator.coordinateVec();
        vec k = integrator.spatialFreqVec();
        vec t = integrator.timeVec();
        vec w = integrator.temporalFreqVec();


        double a = 2.1;
        double wd = w(w.n_elem/2+1);

        OscillatingGaussian S(&integrator, a, wd);

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        S.computeSpatiotemporal();
        g.set_real(S.spatioTemporal());


        //fourier signal
        S.computeFourierTransform();
        f = S.fourierTransform();

        // Backward
        G = integrator.integrate(f);
        G = FFTHelper::fftShift(G);

        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);

                }
            }
        }
    }



    TEST(gaussSpatialCosineTemporal){
        //Mesh
        int ns = 7;
        int nt = 4;
        double ds = 0.1;
        double dt = 0.1;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        Integrator integrator(nt, dt, ns, ds);

        vec s = integrator.coordinateVec();
        vec k = integrator.spatialFreqVec();
        vec t = integrator.timeVec();
        vec w = integrator.temporalFreqVec();

        double a = 2.1;
        double wd = w(w.n_elem/2+3);

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = gauss(a, {s[i], s[j]})
                            * cos(wd * t[l]);
                }
            }
        }


        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = gaussFT(a, {k[i], k[j]})
                            * Functions::delta(wd, w[l])*2*PI;
                }
            }
        }

        f /= integrator.temporalFreqResolution();

        // Backward
        G = integrator.integrate(f);
        G = FFTHelper::fftShift(G);

        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);

                }
            }
        }
    }



    TEST(gaussSpatial){
        //Mesh
        int ns = 7;
        double ds = 0.1;
        double a = 2.1;

        Integrator integrator(0, 0, ns, ds);

        int Ns = pow(2,ns);
        cx_mat g = zeros<cx_mat>(Ns, Ns);
        cx_mat G = zeros<cx_mat>(Ns, Ns);
        cx_mat f = zeros<cx_mat>(Ns, Ns);

        vec s = integrator.coordinateVec();
        vec k = integrator.spatialFreqVec();


        //signal
        for(int i = 0; i < Ns; i++){
            for(int j = 0; j < Ns; j++){
                g(i,j) = gauss(a, {s[i], s[j]});
            }
        }


        //fourier signal
        for(int i = 0; i < Ns; i++){
            for(int j = 0; j < Ns; j++){
                f(i,j) = gaussFT(a, {k[i], k[j]});
            }
        }

        // Backward
        G = integrator.integrate(f);
        G = FFTHelper::fftShift(G);


        //Test
        for(int i = 0; i < Ns; i++){
            for(int j = 0; j < Ns; j++){
                CHECK_CLOSE(real(g(i,j)),
                            real(G(i,j)), 1e-9);

            }
        }
    }



}


