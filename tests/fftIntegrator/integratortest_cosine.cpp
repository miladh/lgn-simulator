#include <unittest++/UnitTest++.h>
#include <armadillo>
#include <iostream>

#include "integrator.h"
#include "stimuli/grating.h"

using namespace std;
using namespace arma;



SUITE(INTEGRATOR){

    TEST(Cosine){
        //Mesh
        int ns = 5;
        int nt = 5;
        double ds = 0.2;
        double dt = 0.1;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        IntegratorSettings settings(nt, dt, ns, ds);
        Integrator integrator(&settings);

        vec s = integrator.coordinateVec();
        vec k = integrator.spatialFreqVec();
        vec t = integrator.timeVec();
        vec w = integrator.temporalFreqVec();


        double wd = w(w.n_elem/2+2);
        double kx = k(k.n_elem/2+6);
        double ky = k(k.n_elem/2+1);


        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = cos(kx*s[i] + ky*s[j] - wd * t[l]);
                }
            }
        }


        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], kx)
                            * Functions::delta(k[j], ky)
                            * Functions::delta(w[l], -wd);
                }
            }
        }

        f *=  8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

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

}



