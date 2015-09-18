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
TEST(fft_ifft){
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
    rowvec f = linspace<rowvec>(-N_2, N_2-1, N) * 1./maxT;
//    f*=2*PI;

    cx_rowvec fSpatial = zeros<cx_rowvec>(N);
    cx_rowvec fSpatial_fftw = zeros<cx_rowvec>(N);
    cx_rowvec fFreq = zeros<cx_rowvec>(N);
//    cx_rowvec fFreq_fftw = zeros<cx_rowvec>(N);


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

    for(int i = 0; i < N; i++){
        fSpatial_fftw(i) *=pow(-1,i) ;
    }


    cout << "----------------------------------------------" << endl;
    cout << "fourier signal: " << endl << abs(fFreq) << endl;
    cout << "----------------------------------------------" << endl;
    cout << "ifft signal: " << endl << real(fSpatial_fftw)<< endl;
    cout << "signal: " << endl << real(fSpatial) << endl;
}


}






SUITE(FFT){

TEST(fft_fb){
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
    cout << "fft signal: " << endl << abs((fFreq_fftw).subvec(0,N_2)) << endl;
    cout << "----------------------------------------------" << endl;
    cout << "ifft signal: " << endl << abs((fSpatial_fftw))*df<< endl;
    cout << "signal: " << endl << real(fSpatial) << endl;
}


}



//TEST(fft){

//    //Spatial Mesh
//    int N = 16;
//    vec mesh = linspace(-0.5, 0.5, N);
//    double dr = mesh(1) - mesh(0);
//    double N_2 = ceil(N/2.);
//    double df = 1./dr/N;
//    double fs = 1./dr;

//    vec f = linspace(-N_2*df, (N - 1. - N_2)*df, N);
//    mesh*=2*PI;

//    DOG dog = DOG(1.0, 1.0, 0., 0.83);

//    cx_mat fSpatial = zeros<cx_mat>(N, N);
//    cx_mat fSpatial_fftw = zeros<cx_mat>(N, N);


//    cx_mat fFreq = zeros<cx_mat>(N, N);
//    cx_mat fFreq_fftw = zeros<cx_mat>(N, N);

//    for(int i = 0; i < N; i++){
//        for(int j = 0; j < N; j++){
//            fSpatial(i,j) = dog.real({mesh[i], mesh[j]});
//            fFreq(i,j) = dog.complex({mesh[i], mesh[j]});
//        }
//    }



//    // Forward
//    int size[2] = {N, N};
//    fftw_complex* in1 = reinterpret_cast<fftw_complex*> (fSpatial.memptr());
//    fftw_complex* out1 = reinterpret_cast<fftw_complex*> (fFreq_fftw.memptr());
//    fftw_plan plan1 = fftw_plan_dft(2, size, in1, out1, FFTW_FORWARD, FFTW_ESTIMATE);

//    fftw_execute(plan1);
//    fftw_destroy_plan(plan1);

//    // Backward
//    fftw_complex* in = reinterpret_cast<fftw_complex*> (fFreq.memptr());
//    fftw_complex* out = reinterpret_cast<fftw_complex*> (fSpatial_fftw.memptr());
//    fftw_plan plan = fftw_plan_dft(2, size, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

//    fftw_execute(plan);
//    fftw_destroy_plan(plan);


//    cout << sum(sum(real(fFreq))) << endl;
//    cout << sum(sum(real(fFreq_fftw)))<< endl;

////    cout << sum(sum(fSpatial)) << endl;
////    cout << sum(sum(real(fSpatial_fftw))) << endl;
////    cout << real(fFreq) << endl;



//}


