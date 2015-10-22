#include <unittest++/UnitTest++.h>
#include <armadillo>
#include <iostream>

#include "integrator.h"
#include "stimuli/patchgrating.h"

using namespace std;
using namespace arma;


SUITE(INTEGRATOR){
    TEST(patchGratingStimuli){
        //Mesh
        int ns = 7;
        int nt = 2;
        double ds = 10.;
        double dt = 0.1;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        Integrator integrator(nt, dt, ns, ds);

        vec s = integrator.coordinateVec();
        vec k = integrator.spatialFreqVec();
        vec t = integrator.timeVec();
        vec w = integrator.temporalFreqVec();

        double C = 1;
        double d = 0.5;
        double wd = w(w.n_elem/2);
        double kx = k(k.n_elem/2);
        double ky = k(k.n_elem/2);
        PatchGrating S(&integrator, {kx, ky}, wd, C, d);

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
//        for(int l = 0; l < Nt; l++){
//            for(int i = 0; i < Ns; i++){
//                for(int j = 0; j < Ns; j++){
//                    CHECK_CLOSE(real(g(i,j,l)),
//                                real(G(i,j,l)), 1e-1);

//                }
//            }
//        }



    }



}


