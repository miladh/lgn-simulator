#include <unittest++/UnitTest++.h>
#include <armadillo>
#include <iostream>

#include "integrator.h"

using namespace std;
using namespace arma;


double gauss(double a, vec r){
    return exp(-a*dot(r,r));
}

double gaussFT(double a, vec k){
    return PI/a*exp(-dot(k,k)/4./a);
}

SUITE(INTEGRATOR){

    TEST(gaussfft){
        //Mesh
        int ns = 7;
        double ds = 0.1;
        double a = 2.1;

        IntegratorSettings settings(0, 0, ns, ds);
        Integrator integrator(&settings);

        int Ns = pow(2,ns);
        cx_mat g = zeros<cx_mat>(Ns, Ns);
        cx_mat G = zeros<cx_mat>(Ns, Ns);
        cx_mat f = zeros<cx_mat>(Ns, Ns);

        vec s = integrator.coordinateVec();
        vec k = integrator.spatialFreqVec();

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

        for(int i = 0; i < Ns; i++){
            for(int j = 0; j < Ns; j++){
                CHECK_CLOSE(real(g(i,j)),
                            real(G(i,j)), 1e-9);

            }
        }
    }



}


