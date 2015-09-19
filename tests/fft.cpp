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
TEST(ifft_of_analytic_function){
    cout << "fft ifft test" << endl;
    cout << "-------------" << endl;
    double pi = acos(-1);

    //Spatial Mesh
    int N = 16;
    double maxT = 1.;
    double g = 2;

    double dt = maxT/N;
    double df = 1./maxT;
    double fs = N/maxT;

    double N_2 = ceil(N/2.);

    rowvec t = linspace<rowvec>(0, maxT-dt, N);
    rowvec f1 = linspace<rowvec>(0, N_2-1, N_2);
    rowvec f2 = linspace<rowvec>(-N_2,-f1[1], N_2);
    rowvec f = join_rows(f1,f2)* 1./maxT;

    cx_rowvec fSpatial = zeros<cx_rowvec>(N);
    cx_rowvec fSpatial_fftw = zeros<cx_rowvec>(N);
    cx_rowvec fFreq = zeros<cx_rowvec>(N);


    cout << "N_2: " << N_2 << endl;
    cout << "dt: " << dt << endl;
    cout << "df: " << df << endl;
    cout << "sampling freq: " << fs << endl;
    cout << "signal freq: "   << g << endl;
    cout << "----------------------------------------------" << endl;

    cout << "Temporal:" << endl << t << endl;
    cout << "Freq:" << endl << f << endl;



    //signal
    for(int i = 0; i < N; i++){
        fSpatial(i) = 0.0 + cos(2*pi*t[i] * g);
    }

    //fourier signal
    for(int i = 0; i < N; i++){
        fFreq(i) = 0.5*(Functions::delta(f[i],g) + Functions::delta(f[i], -g));
    }


//    fFreq *= df;


    // Backward
    int size[1] = {N};
    fftw_complex* in = reinterpret_cast<fftw_complex*> (fFreq.memptr());
    fftw_complex* out = reinterpret_cast<fftw_complex*> (fSpatial_fftw.memptr());
    fftw_plan plan = fftw_plan_dft(1, size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_execute(plan);
    fftw_destroy_plan(plan);

    fSpatial_fftw *=df;

    cout << "----------------------------------------------" << endl;
    cout << "fourier signal: " << endl << real(fFreq) << endl;
    cout << "----------------------------------------------" << endl;
    cout << "ifft signal: " << endl << real(fSpatial_fftw)<< endl;
    cout << "signal: " << endl << real(fSpatial) << endl;
}


}






SUITE(FFT){

TEST(fft_of_analytic_function_followed_by_ifft){
    double pi = acos(-1);


    //Spatial Mesh
    int N = 16;
    double maxT = 1.;
    double g = 2;

    double dt = maxT/N;
    double df = 1./maxT;
    double fs = N/maxT;


    double N_2 = ceil(N/2.);

    rowvec t = linspace<rowvec>(0, maxT-dt, N);
    rowvec f = linspace<rowvec>(0, N-1, N) * 1./maxT;
//    f*=2*PI;

    cx_rowvec fSpatial = zeros<cx_rowvec>(N);
    cx_rowvec fSpatial_fftw = zeros<cx_rowvec>(N);
    cx_rowvec fFreq = zeros<cx_rowvec>(N);
    cx_rowvec fFreq_fftw = zeros<cx_rowvec>(N);


    cout << "dt: " << dt << endl;
    cout << "df: " << df << endl;
    cout << "sampling freq: " << fs << endl;
    cout << "signal freq: "   << g << endl;
    cout << "----------------------------------------------" << endl;

    cout << "Temporal:" << endl << t << endl;
    cout << "Freq:" << endl << f.subvec(0,N_2) << endl;


    //signal
    for(int i = 0; i < N; i++){
        fSpatial(i) = 0.0 + cos(2*pi*t[i] * g);
    }


    // Forward
    int size[1] = {N};
    fftw_complex* in1 = reinterpret_cast<fftw_complex*> (fSpatial.memptr());
    fftw_complex* out1 = reinterpret_cast<fftw_complex*> (fFreq_fftw.memptr());
    fftw_plan plan1 = fftw_plan_dft(1, size, in1, out1, FFTW_FORWARD,FFTW_ESTIMATE);

    fftw_execute(plan1);
    fftw_destroy_plan(plan1);

    fFreq_fftw *= dt;

    // Backward
    fftw_complex* in = reinterpret_cast<fftw_complex*> (fFreq_fftw.memptr());
    fftw_complex* out = reinterpret_cast<fftw_complex*> (fSpatial_fftw.memptr());
    fftw_plan plan = fftw_plan_dft(1, size, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(plan);
    fftw_destroy_plan(plan);


    cout << "----------------------------------------------" << endl;
    cout << "fft signal: " << endl << real((fFreq_fftw)) << endl;
    cout << "----------------------------------------------" << endl;
    cout << "ifft signal: " << endl << real((fSpatial_fftw))*df<< endl;
    cout << "signal: " << endl << real(fSpatial) << endl;
}


}



