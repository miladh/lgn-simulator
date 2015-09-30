#include <unittest++/UnitTest++.h>
#include <armadillo>
#include <iostream>


<<<<<<< HEAD
#include "stimuli/patchgrating.h"
=======
#include "integrator.h"
#include "spatialKernels/dog.h"
>>>>>>> fft_3d


using namespace std;
using namespace arma;


SUITE(INTEGRATOR){

    TEST(dogfft){
        //Mesh
        int nt = 0;
        int ns = 4;
        double maxT = 1.;

        DOG dog(1.0, 10.5, 0.0, 1.0);

        IntegratorSettings settings(nt,ns,maxT);
        Integrator integrator(&settings);


        int Ns = pow(2,ns);
        cx_mat fSpatial = zeros<cx_mat>(Ns, Ns);
        cx_mat fSpatial_fftw = zeros<cx_mat>(Ns, Ns);
        cx_mat fFreq = zeros<cx_mat>(Ns, Ns);


        vec s = integrator.coordinateVec();
        vec k = integrator.spatialFreqVec();

        cout << s.t() << endl;
        cout << k.t() << endl;


        for(int i = 0; i < Ns; i++){
            for(int j = 0; j < Ns; j++){
                fSpatial(i,j) = dog.real({s[i], s[j]});
            }
        }


        //fourier signal
        for(int i = 0; i < Ns; i++){
            for(int j = 0; j < Ns; j++){
                fFreq(i,j) = dog.complex({k[i], k[j]});
            }
        }


        // Backward
        fSpatial_fftw = integrator.integrate(fFreq);
        fSpatial_fftw = FFTHelper::fftShift(fSpatial_fftw);

//        for(int i = 0; i < Ns; i++){
//            for(int j = 0; j < Ns; j++){
//                CHECK_CLOSE(real(fSpatial(i,j)),
//                            real(fSpatial_fftw(i,j)), 1e-8);

//            }
//        }

        cout << real(fFreq) << endl;
        cout << real(fSpatial) << endl;
        cout << real(fSpatial_fftw) << endl;


    }



}


