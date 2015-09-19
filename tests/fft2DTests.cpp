#include <unittest++/UnitTest++.h>
#include <armadillo>
#include <iostream>
#include <fftw3.h>

#include "stimuli/stimuli.h"
#include "spatialKernels/dog.h"
#include "math/functions.h"


using namespace std;
using namespace arma;




SUITE(DEVELOPMENT){

    /************************************************
     * FFTW_BACKWARD of 2D analytic function
     * */
    TEST(ifft){
        cout <<endl << "ifft 2d" << endl;
        cout << "----------------------" << endl;
        double pi = acos(-1);

        //Spatial Mesh
        int N = 8;
        double maxT = 1.;
        double gx = 1;
        double gy = 2;

        double dt = maxT/N;
        double df = 1./maxT;
        double fs = N/maxT;

        double N_2 = ceil(N/2.);

        rowvec t = linspace<rowvec>(0, maxT-dt, N);
        rowvec f1 = linspace<rowvec>(0, N_2-1, N_2);
        rowvec f2 = linspace<rowvec>(-N_2,-f1[1], N_2);
        rowvec f = join_rows(f1,f2)* 1./maxT;

        cx_mat fSpatial = zeros<cx_mat>(N, N);
        cx_mat fSpatial_fftw = zeros<cx_mat>(N, N);
        cx_mat fFreq = zeros<cx_mat>(N, N);


        cout << "N_2: " << N_2 << endl;
        cout << "dt: " << dt << endl;
        cout << "df: " << df << endl;
        cout << "sampling freq: " << fs << endl;
        cout << "signal freq: "   << gx << " - " << gy << endl;
        cout << "----------------------------------------------" << endl;

        cout << "Temporal:" << endl << t << endl;
        cout << "Freq:" << endl << f << endl;


        //signal
        for(int i = 0; i < N; i++){
            for(int j = 0; j < N; j++){
                fSpatial(i,j) = cos( 2*pi* (t[i] * gx + t[j] * gy) );
            }
        }

        //fourier signal
        for(int i = 0; i < N; i++){
            for(int j = 0; j < N; j++){
                fFreq(i,j) = 0.5*(
                            Functions::delta(f[i],gx)* Functions::delta(f[j],gy)
                            +Functions::delta(f[i], -gx)*Functions::delta(f[j], -gy)
                            );
            }
        }


        // Backward
        int size[2] = {N, N};
        fftw_complex* in = reinterpret_cast<fftw_complex*> (fFreq.memptr());
        fftw_complex* out = reinterpret_cast<fftw_complex*> (fSpatial_fftw.memptr());
        fftw_plan plan = fftw_plan_dft(2, size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

        fftw_execute(plan);
        fftw_destroy_plan(plan);

        fSpatial_fftw *=df*df;

        cout << "----------------------------------------------" << endl;
        cout << "fourier signal: " << endl << real(fFreq) << endl;
        cout << "----------------------------------------------" << endl;
        cout << "ifft signal: " << endl << real(fSpatial_fftw)<< endl;
        cout << "signal: " << endl << real(fSpatial) << endl;


        for(int i = 0; i < N; i++){
            for(int j = 0; j < N; j++){
                CHECK_CLOSE(real(fSpatial_fftw(i,j)), real(fSpatial(i,j)), 1e-8);

            }
        }
    }


}



