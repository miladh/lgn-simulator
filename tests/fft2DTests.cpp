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
     * FFTW_BACKWARD of 3D analytic function
     * */
    TEST(ifft_3d){
        cout <<endl << "ifft 3d" << endl;
        cout << "----------------------" << endl;
        double pi = acos(-1);

        //Temporal Mesh
        int Nt = 8;
        double maxT = 1.;
        double wd = 2;
        double dt = maxT/Nt;
        double dw = 1./maxT;
        double ws = Nt/maxT;
        double Nt_2 = ceil(Nt/2.);

        rowvec t = linspace<rowvec>(0, maxT-dt, Nt);
        rowvec w1 = linspace<rowvec>(0, Nt_2-1, Nt_2);
        rowvec w2 = linspace<rowvec>(-Nt_2,-w1[1], Nt_2);
        rowvec w = join_rows(w1,w2)* 1./maxT;

        //Spatial Mesh x
        int Nx = 8;
        double maxX = 1.;
        double kdx = 3;
        double dx = maxX/Nx;
        double dkx = 1./maxX;
        double kxs = Nx/maxX;
        double Nx_2 = ceil(Nx/2.);

        rowvec x = linspace<rowvec>(0, maxX-dx, Nx);
        rowvec kx1 = linspace<rowvec>(0, Nx_2-1, Nx_2);
        rowvec kx2 = linspace<rowvec>(-Nx_2,-kx1[1], Nx_2);
        rowvec kx = join_rows(kx1,kx2)* 1./maxX;


        //Spatial Mesh y
        int Ny = 4;
        double maxY = 1.;
        double kdy = 1;
        double dy = maxY/Ny;
        double dky = 1./maxY;
        double kys = Nx/maxY;
        double Ny_2 = ceil(Ny/2.);

        rowvec y = linspace<rowvec>(0, maxY-dy, Ny);
        rowvec ky1 = linspace<rowvec>(0, Ny_2-1, Ny_2);
        rowvec ky2 = linspace<rowvec>(-Ny_2,-ky1[1], Ny_2);
        rowvec ky = join_rows(ky1,ky2)* 1./maxY;

        cx_cube fSpatial = zeros<cx_cube>(Nx, Ny, Nt);
        cx_cube fSpatial_fftw = zeros<cx_cube>(Nx, Ny, Nt);
        cx_cube fFreq = zeros<cx_cube>(Nx, Ny, Nt);

        cout << "dt: " << dt << ", dx: " << dx << ", dy: " << dy << endl;
        cout << "dw: " << dw << ", dkx: " << dkx  << ", dky: " << dky << endl;
        cout << "sampling ws: " << ws  <<
                ", sampling kxs: " << kxs <<
                ", sampling kys: " << kys<< endl;
        cout << "signal wd: "   << wd
             << ", signal kdx: " << kdx
             << ", signal kdy: " << kdy << endl;
        cout << "----------------------------------------------" << endl;

        cout << "Temporal:" << endl << t << endl;
        cout << "Spatial x:" << endl << x << endl;
        cout << "Spatial y:" << endl << y << endl;
        cout << "----------------------------------------------" << endl;
        cout << "Freq temporal:" << endl << w << endl;
        cout << "Freq spatial x:" << endl << kx << endl;
        cout << "Freq spatial y:" << endl << ky << endl;


        //signal
        for(int k = 0; k < Nt; k++){
            for(int i = 0; i < Nx; i++){
                for(int j = 0; j < Ny; j++){
                    fSpatial(i,j, k) = cos( 2*pi* (kdx * x[i]+ kdy * y[j]
                                                   - wd * t[k]) );
                }
            }
        }

        //fourier signal
        for(int k = 0; k < Nt; k++){
            for(int i = 0; i < Nx; i++){
                for(int j = 0; j < Ny; j++){
                    fFreq(i,j,k) = 0.5*(
                                Functions::delta(kx[i],kdx)
                                *Functions::delta(ky[j],kdy)
                                *Functions::delta(w[k],-wd)
                                +Functions::delta(kx[i],-kdx)
                                *Functions::delta(ky[j],-kdy)
                                *Functions::delta(w[k],wd)
                                );
                }
            }
        }

        // Backward
        int size[3] = {Nx, Ny, Nt};
        fftw_complex* in = reinterpret_cast<fftw_complex*> (fFreq.memptr());
        fftw_complex* out = reinterpret_cast<fftw_complex*> (fSpatial_fftw.memptr());
        fftw_plan plan = fftw_plan_dft(3, size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

        fftw_execute(plan);
        fftw_destroy_plan(plan);

        fSpatial_fftw *=dkx*dky*dw;


//        cout << "----------------------------------------------" << endl;
//        cout << "fourier signal: " << endl << real(fFreq) << endl;
//        cout << "----------------------------------------------" << endl;
//        cout << "ifft signal: " << endl << real(fSpatial_fftw)<< endl;
//        cout << "signal: " << endl << real(fSpatial) << endl;


        for(int k = 0; k < Nt; k++){
            for(int i = 0; i < Nx; i++){
                for(int j = 0; j < Ny; j++){
                    CHECK_CLOSE(real(fSpatial(i,j,k)),
                                real(fSpatial_fftw(i,j,k)), 1e-8);

                }
            }
        }


    }


}


/************************************************
 * FFTW_BACKWARD of 2D analytic function
 * */
TEST(ifft_2d){
    cout <<endl << "ifft 2d" << endl;
    cout << "----------------------" << endl;
    double pi = acos(-1);

    //Temporal Mesh
    int Nt = 8;
    double maxT = 1.;
    double wd = 2;
    double dt = maxT/Nt;
    double dw = 1./maxT;
    double ws = Nt/maxT;
    double Nt_2 = ceil(Nt/2.);

    rowvec t = linspace<rowvec>(0, maxT-dt, Nt);
    rowvec w1 = linspace<rowvec>(0, Nt_2-1, Nt_2);
    rowvec w2 = linspace<rowvec>(-Nt_2,-w1[1], Nt_2);
    rowvec w = join_rows(w1,w2)* 1./maxT;

    //Spatial Mesh
    int Nx = 4;
    double maxX = 1.;
    double kd = 1;
    double dx = maxX/Nx;
    double dk = 1./maxX;
    double ks = Nx/maxX;
    double Nx_2 = ceil(Nx/2.);

    rowvec x = linspace<rowvec>(0, maxX-dx, Nx);
    rowvec k1 = linspace<rowvec>(0, Nx_2-1, Nx_2);
    rowvec k2 = linspace<rowvec>(-Nx_2,-k1[1], Nx_2);
    rowvec k = join_rows(k1,k2)* 1./maxX;


    cx_mat fSpatial = zeros<cx_mat>(Nt, Nx);
    cx_mat fSpatial_fftw = zeros<cx_mat>(Nt, Nx);
    cx_mat fFreq = zeros<cx_mat>(Nt, Nx);

    cout << "dt: " << dt << ", dx: " << dx << endl;
    cout << "dw: " << dw << ", dk: " << dk  << endl;
    cout << "sampling ws: " << ws  << ", sampling ks: " << ks << endl;
    cout << "signal wd: "   << wd << ", signal kd: " << kd << endl;
    cout << "----------------------------------------------" << endl;

    cout << "Temporal:" << endl << t << endl;
    cout << "Spatial:" << endl << x << endl;
    cout << "----------------------------------------------" << endl;
    cout << "Freq temporal:" << endl << w << endl;
    cout << "Freq spatial:" << endl << k << endl;


    //signal
    for(int i = 0; i < Nt; i++){
        for(int j = 0; j < Nx; j++){
            fSpatial(i,j) = cos( 2*pi* (kd * x[j] - wd * t[i]) );
        }
    }

    //fourier signal
    for(int i = 0; i < Nt; i++){
        for(int j = 0; j < Nx; j++){
            fFreq(i,j) = 0.5*(
                        Functions::delta(k[j],kd)* Functions::delta(w[i],-wd)
                        +Functions::delta(k[j], -kd)*Functions::delta(w[i], wd)
                        );
        }
    }

    // Backward
    int size[2] = {Nx, Nt};
    fftw_complex* in = reinterpret_cast<fftw_complex*> (fFreq.memptr());
    fftw_complex* out = reinterpret_cast<fftw_complex*> (fSpatial_fftw.memptr());
    fftw_plan plan = fftw_plan_dft(2, size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_execute(plan);
    fftw_destroy_plan(plan);

    fSpatial_fftw *=dk*dw;


    cout << "----------------------------------------------" << endl;
    cout << "fourier signal: " << endl << real(fFreq) << endl;
    cout << "----------------------------------------------" << endl;
    cout << "ifft signal: " << endl << real(fSpatial_fftw)<< endl;
    cout << "signal: " << endl << real(fSpatial) << endl;


    for(int i = 0; i < Nt; i++){
        for(int j = 0; j < Nx; j++){
            CHECK_CLOSE(real(fSpatial(i,j)), real(fSpatial_fftw(i,j)), 1e-8);

        }
    }
}



