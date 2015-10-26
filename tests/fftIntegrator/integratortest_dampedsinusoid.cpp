//#include <unittest++/UnitTest++.h>
//#include <armadillo>
//#include <iostream>

//#include "integrator.h"
//#include "math/functions.h"

//#include "stimuli/oscillatinggaussian.h"

//using namespace std;
//using namespace arma;


//double dampedSinusoid(double a, double t0, double wd, vec r, double t){

//    double td = t - t0;
//    return exp(-a*td) /** sin(wd * td)*/;
//}

//complex<double> dampedSinusoidFT(double a, double wd, vec k, double w){
//    complex<double> i = {0,1};
////    complex<double> g = wd / (-w*w + 2*w*a*i + wd*wd + a*a);
//    complex<double> g = 2.*a / (a*a + w*w);

//    return g;
//}

//SUITE(INTEGRATOR){

//    TEST(dampedSinusoid){
//        //Mesh
//        int ns = 8;
//        int nt = 8;
//        double ds = 0.1;
//        double dt = 0.1;

//        int Ns = pow(2,ns);
//        int Nt = pow(2,nt);

//        Integrator integrator(nt, dt, ns, ds);

//        vec s = integrator.coordinateVec();
//        vec k = integrator.spatialFreqVec();
//        vec t = integrator.timeVec();
//        vec w = integrator.temporalFreqVec();

//        double a = 10.0;
//        double wd = w(w.n_elem/2+1);
//        double t0 = t.min();


//        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
//        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
//        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


//        //Spatiotemporal signal
//        for(int l = 0; l < Nt; l++){
//            for(int i = 0; i < Ns; i++){
//                for(int j = 0; j < Ns; j++){
//                    g(i,j,l) = dampedSinusoid(a,t0, wd, {s[i], s[j]}, t[l]);
//                }
//            }
//        }


//        //fourier signal
//        for(int l = 0; l < Nt; l++){
//            for(int i = 0; i < Ns; i++){
//                for(int j = 0; j < Ns; j++){
//                    f(i,j,l) = dampedSinusoidFT(a, wd, {k[i], k[j]}, w[l]);
//                    f(i,j,l) *= Functions::delta(k[i],0)*Functions::delta(k[j],0)
//                            * 4*PI*PI;
//                }
//            }
//        }

//        f /= integrator.spatialFreqResolution()/integrator.spatialFreqResolution();

//        // Backward
//        G = integrator.integrate(f);
//        G = FFTHelper::fftShift(G);

//        // Test
////        for(int l = 0; l < Nt; l++){
////            for(int i = 0; i < Ns; i++){
////                for(int j = 0; j < Ns; j++){
////                    CHECK_CLOSE(real(g(i,j,l)),
////                                real(G(i,j,l)), 1e-9);

////                }
////            }
////        }
//    }



//}



