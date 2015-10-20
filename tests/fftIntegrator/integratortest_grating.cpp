//#include <unittest++/UnitTest++.h>
//#include <armadillo>
//#include <iostream>

//#include "integrator.h"
//#include "stimuli/grating.h"

//using namespace std;
//using namespace arma;



//SUITE(INTEGRATOR){
//    TEST(gratingStimuli){
//        //Mesh
//        int ns = 9;
//        int nt = 1;
//        double ds = 0.1;
//        double dt = 0.1;

//        int Ns = pow(2,ns);
//        int Nt = pow(2,nt);


//        IntegratorSettings settings(nt, dt, ns, ds);
//        Integrator integrator(&settings);

//        vec s = integrator.coordinateVec();
//        vec k = integrator.spatialFreqVec();
//        vec t = integrator.timeVec();
//        vec w = integrator.temporalFreqVec();

//        double C = 1.0;
//        double wd = w[w.n_elem/2+1];
//        vec kd = {k[k.n_elem/2+3], 0.0};
//        Grating S(integrator, kd, wd, C);

//        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
//        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
//        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


//        //Spatiotemporal signal
//        S.computeSpatiotemporal();
//        g.set_real(S.spatioTemporal());


//        //fourier signal
//        //Spatiotemporal signal
//        S.computeFourierTransform();
//        f = S.fourierTransform();


//        // Backward
//        G = integrator.integrate(f);
//        G = FFTHelper::fftShift(G);

//        // Test
//        for(int l = 0; l < Nt; l++){
//            for(int i = 0; i < Ns; i++){
//                for(int j = 0; j < Ns; j++){
//                    CHECK_CLOSE(real(g(i,j,l)),
//                                real(G(i,j,l)), 1e-9);

//                }
//            }
//        }



//    }



//}


