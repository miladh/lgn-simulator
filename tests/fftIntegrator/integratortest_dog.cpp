#include <unittest++/UnitTest++.h>
#include <armadillo>
#include <iostream>

#include "integrator.h"
#include "spatialKernels/dog.h"


using namespace std;
using namespace arma;


SUITE(INTEGRATOR){

    TEST(dogfft){
        //Mesh
        int ns = 8;
        double ds = 0.1;
        double A = 1.0;
        double a = 2.1;

        IntegratorSettings settings(0, 0, ns, ds);
        Integrator integrator(&settings);

        int Ns = pow(2,ns);
        cx_mat g = zeros<cx_mat>(Ns, Ns);
        cx_mat G = zeros<cx_mat>(Ns, Ns);
        cx_mat f = zeros<cx_mat>(Ns, Ns);

        vec s = integrator.coordinateVec();
        vec k = integrator.spatialFreqVec();

        DOG dog(A, a, 0, 0.1);

        for(int i = 0; i < Ns; i++){
            for(int j = 0; j < Ns; j++){
                g(i,j) = dog.spatial({s[i], s[j]});
            }
        }


        //fourier signal
        for(int i = 0; i < Ns; i++){
            for(int j = 0; j < Ns; j++){
                f(i,j) = dog.fourierTransform({k[i], k[j]});
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


