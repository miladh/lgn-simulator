#include <unittest++/UnitTest++.h>
#include <armadillo>
#include <iostream>
#include <fftw3.h>

#include "stimuli/stimuli.h"
#include "spatialKernels/dog.h"
#include "math/functions.h"


using namespace std;
using namespace arma;




SUITE(FFT_1D){

    /************************************************
     * FFTW_FORWARD of analytic function, multiplied by
     * another analytic function, followed
     * by FFTW_BACKWARD of the product.
     * */
    TEST(fft_multiply_ifft){
        cout <<endl << "fft - multiply - ifft" << endl;
        cout << "----------------------" << endl;
        double pi = acos(-1);

        //Spatial Mesh
        int N = 16;
        double maxT = 1.;

        // Signals:
        double g1 = 2;
        double g2 = 4;

        double dt = maxT/N;
        double df = 1./maxT;
        double fs = N/maxT;

        double N_2 = ceil(N/2.);

        rowvec t = linspace<rowvec>(0, maxT-dt, N);
        rowvec f1 = linspace<rowvec>(0, N_2-1, N_2);
        rowvec f2 = linspace<rowvec>(-N_2,-f1[1], N_2);
        rowvec f = join_rows(f1,f2)* 1./maxT;


        cx_rowvec signal1 = zeros<cx_rowvec>(N);
        cx_rowvec signal1_ft = zeros<cx_rowvec>(N);
        cx_rowvec signal2_ft = zeros<cx_rowvec>(N);

        cx_rowvec fSpatial = zeros<cx_rowvec>(N);
        cx_rowvec fFreq = zeros<cx_rowvec>(N);
        cx_rowvec fSpatial_fftw = zeros<cx_rowvec>(N);


        cout << "N_2: " << N_2 << endl;
        cout << "dt: " << dt << endl;
        cout << "df: " << df << endl;
        cout << "sampling freq: " << fs << endl;
        cout << "signal 1 freq: "   << g1 << endl;
        cout << "signal 2 freq: "   << g2 << endl;
        cout << "----------------------------------------------" << endl;

        cout << "Temporal:" << endl << t << endl;
        cout << "Freq:" << endl << f << endl;


        // signal1 + signal2
        // signal1
        // fourier signal2
        for(int i = 0; i < N; i++){
            fSpatial(i) = cos(2*pi*t[i] * g1) + cos(2*pi*t[i] * g2);
            signal1(i) = cos(2*pi*t[i] * g1);
            signal2_ft(i)=0.5*(Functions::delta(f[i],g2)+Functions::delta(f[i],-g2));
        }

        // Forward of signal 1
        int size[1] = {N};
        fftw_complex* in1 = reinterpret_cast<fftw_complex*> (signal1.memptr());
        fftw_complex* out1 = reinterpret_cast<fftw_complex*> (signal1_ft.memptr());
        fftw_plan plan1=fftw_plan_dft(1, size, in1, out1, FFTW_FORWARD,FFTW_ESTIMATE);
        fftw_execute(plan1);
        fftw_destroy_plan(plan1);

        signal1_ft *= dt;

        //adding the fourier signals:
        fFreq = signal1_ft + signal2_ft;

        // Backward of product of signal 1 and signal 2
        fftw_complex* in = reinterpret_cast<fftw_complex*> (fFreq.memptr());
        fftw_complex* out = reinterpret_cast<fftw_complex*> (fSpatial_fftw.memptr());
        fftw_plan plan = fftw_plan_dft(1, size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);

        fSpatial_fftw *=df;

        cout << "----------------------------------------------" << endl;
        cout << "fourier signal1: " << endl << real(signal1_ft) << endl;
        cout << "fourier signal2: " << endl << real(signal2_ft) << endl;
        cout << "fourier signal: " << endl << real(fFreq) << endl;
        cout << "----------------------------------------------" << endl;
        cout << "ifft signal: " << endl << real(fSpatial_fftw)<< endl;
        cout << "signal: " << endl << real(fSpatial) << endl;


        for(int i = 0; i < N; i++){
            CHECK_CLOSE(real(fSpatial_fftw(i)), real(fSpatial(i)), 1e-8);
        }
    }


    /************************************************
     * FFTW_BACKWARD of analytic function
     * */
    TEST(ifft){
        cout <<endl << "ifft" << endl;
        cout << "----------------------" << endl;
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

        for(int i = 0; i < N; i++){
            CHECK_CLOSE(real(fSpatial_fftw(i)), real(fSpatial(i)), 1e-8);
        }
    }


    /************************************************
     * FFTW_FORWARD of analytic function, followed
     * by FFTW_BACKWARD.
     * */
    TEST(fft_ifft){
        cout <<endl<< "fft - ifft" << endl;
        cout << "----------------------" << endl;
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

        for(int i = 0; i < N; i++){
            CHECK_CLOSE(real(fSpatial_fftw(i)), real(fSpatial(i)), 1e-8);
        }
    }




}



