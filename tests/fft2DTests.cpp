#include <unittest++/UnitTest++.h>
#include <armadillo>
#include <iostream>
#include <fftw3.h>

#include "math/functions.h"
#include "integrator.h"


using namespace std;
using namespace arma;


double gaussian(vec f, double x, double y, double t){

    double wd = f[0];
    double kdx =f[1];
    double kdy = f[2];
    double g = cos(2*PI*(kdx * x+ kdy * y - wd * t) );

    return g;
}


double fft_gaussian(vec f, double kx, double ky, double w){
    double wd = f[0];
    double kdx =f[1];
    double kdy = f[2];

    double g =  0.5*(
                Functions::delta(kx,kdx*2*PI)
                *Functions::delta(ky,kdy*2*PI)
                *Functions::delta(w,-wd*2*PI)
                +Functions::delta(kx,-kdx*2*PI)
                *Functions::delta(ky,-kdy*2*PI)
                *Functions::delta(w,wd*2*PI)
                );

    return g;
}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

SUITE(FFT_nD){

    /************************************************
     * FFTW_BACKWARD of 3D analytic function
     * */
    TEST(ifft_with_integrator){
        double pi = acos(-1);

        //Mesh
        int nt = 5;
        int ns = 3;
        double maxT = 2.;

        IntegratorSettings settings(nt,ns,maxT);
        Integrator integrator(&settings);

        int Nt = pow(2,nt);
        int Ns = pow(2,ns);
        cx_cube fSpatial = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube fSpatial_fftw = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube fFreq = zeros<cx_cube>(Ns, Ns, Nt);

        vec t = integrator.timeVec();
        vec s = integrator.coordinateVec();
        vec w = integrator.temporalFreqVec(); // FACTOR 2PI!!!!
        vec k = integrator.spatialFreqVec(); // FACTOR 2PI!!!!

        //signal
        vec f = {2,1,1};
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    fSpatial(i,j, l) = gaussian(f, s[i], s[j], t[l]);
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    fFreq(i,j,l) = fft_gaussian(f, k[i], k[j], w[l]);
                }
            }
        }

        // Backward
        fSpatial_fftw = integrator.integrate(fFreq);

        for(int k = 0; k < Nt; k++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(fSpatial(i,j,k)),
                                real(fSpatial_fftw(i,j,k)), 1e-8);

                }
            }
        }


//        cout << "Temporal:" << endl << t.t() << endl;
//        cout << "Spatial:" << endl << s.t() << endl;
//        cout << "----------------------------------------------" << endl;
//        cout << "Freq temporal:" << endl << w.t() << endl;
//        cout << "Freq spatial:" << endl << k.t() << endl;
//        cout << "----------------------------------------------" << endl;
//        cout << "fourier signal: " << endl << real(fFreq) << endl;
//        cout << "----------------------------------------------" << endl;
//        cout << "ifft signal: " << endl << real(fSpatial_fftw)<< endl;
//        cout << "signal: " << endl << real(fSpatial) << endl;


    }



}




