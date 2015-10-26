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
                                real(G(i,j,l)), 1e-12);

                }
            }
        }
    }



    TEST(gauss_0){
        int ns = 9;
        int nt = 2;
        double ds = 0.1;
        double dt = 0.1;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        Integrator integrator(nt, dt, ns, ds);

        vec s = integrator.coordinateVec();
        vec k = integrator.spatialFreqVec();
        vec t = integrator.timeVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);

        double wd = w(2);
        double a =0.1;

        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) =gauss(a, {s[i], s[j]}) * cos(wd * t[l]);
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = gaussFT(a, {k[i], k[j]})
                            * Functions::delta(wd, w[l]);
                }
            }
        }

        f *= 2*PI;
        f /= integrator.temporalFreqResolution();

        // Backward
        G = integrator.integrate(f);
        G = FFTHelper::fftShift(G);

        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-12);
                }
            }
        }
    }
    TEST(gauss_1){
        int ns = 8;
        int nt = 2;
        double ds = 0.1;
        double dt = 0.1;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        Integrator integrator(nt, dt, ns, ds);

        vec s = integrator.coordinateVec();
        vec k = integrator.spatialFreqVec();
        vec t = integrator.timeVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);

        double wd = w(2);
        double a =2.2;

        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) =gauss(a, {s[i], s[j]}) * cos(wd * t[l]);
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = gaussFT(a, {k[i], k[j]})
                            * Functions::delta(wd, w[l]);
                }
            }
        }

        f *= 2*PI;
        f /= integrator.temporalFreqResolution();

        // Backward
        G = integrator.integrate(f);
        G = FFTHelper::fftShift(G);

        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-12);
                }
            }
        }
    }
    TEST(gauss_2){
        int ns = 8;
        int nt = 2;
        double ds = 0.1;
        double dt = 0.1;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        Integrator integrator(nt, dt, ns, ds);

        vec s = integrator.coordinateVec();
        vec k = integrator.spatialFreqVec();
        vec t = integrator.timeVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);

        double wd = w(2);
        double a =5.8;

        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) =gauss(a, {s[i], s[j]}) * cos(wd * t[l]);
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = gaussFT(a, {k[i], k[j]})
                            * Functions::delta(wd, w[l]);
                }
            }
        }

        f *= 2*PI;
        f /= integrator.temporalFreqResolution();

        // Backward
        G = integrator.integrate(f);
        G = FFTHelper::fftShift(G);

        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-12);
                }
            }
        }
    }
    TEST(gauss_3){
        int ns = 8;
        int nt = 2;
        double ds = 0.01;
        double dt = 0.1;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        Integrator integrator(nt, dt, ns, ds);

        vec s = integrator.coordinateVec();
        vec k = integrator.spatialFreqVec();
        vec t = integrator.timeVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);

        double wd = w(2);
        double a =60.8;

        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) =gauss(a, {s[i], s[j]}) * cos(wd * t[l]);
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = gaussFT(a, {k[i], k[j]})
                            * Functions::delta(wd, w[l]);
                }
            }
        }

        f *= 2*PI;
        f /= integrator.temporalFreqResolution();

        // Backward
        G = integrator.integrate(f);
        G = FFTHelper::fftShift(G);

        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-12);
                }
            }
        }
    }


}


