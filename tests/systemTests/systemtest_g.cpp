#include <unittest++/UnitTest++.h>

#include "neurons/ganglioncell.h"
#include "neurons/relaycell.h"

#include "stimuli/grating.h"
#include "integrator.h"
#include "spatialKernels/dog.h"
#include "temporalKernels/temporallyconstant.h"
#include "temporalKernels/decayingexponential.h"



SUITE(SYSTEM){


    TEST(responseDecayingDOG){
        int ns = 5;
        int nt = 5;
        double ds = 0.1;
        double dt = 0.1;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double A = 1.0;
        double a = 2.1;
        double B = 0.5;
        double b = 0.8;

        double tau = 1.0;
        double delay = 0.0;

        double C = -2.3;

        cube R = zeros<cube>(Ns, Ns, Nt);
        cube Rex = zeros<cube>(Ns, Ns, Nt);

        //Integrator
        Integrator integrator(nt, dt, ns, ds);
        vec s = integrator.coordinateVec();
        vec k = integrator.spatialFreqVec();
        vec t = integrator.timeVec();
        vec w = integrator.temporalFreqVec();


        //Cell
        DOG dog(A, a, B, b);
        DecayingExponential Kt(tau, delay);
        GanglionCell ganglion(&integrator, &dog, &Kt);

        //Stimulus
        double wd = w(2);
        double kx = k(1);
        double ky = k(1);
        Grating S(&integrator, {kx, ky}, wd, C);
        S.computeFourierTransform();


        //Compute analytic:
        double W = dog.fourierTransform({kx, ky}) * Kt.fourierTransform(wd);
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    Rex(i,j,l) = C * W * cos(kx*s[i] + ky*s[j] - wd * t[l]);

                }
            }
        }


        //Compute numerical
        ganglion.computeResponse(&S);
        R = ganglion.response();


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(Rex(i,j,l), R(i,j,l), 1e-12);

                }
            }
        }

    }



    TEST(responseDOG){
        int ns = 5;
        int nt = 5;
        double ds = 0.1;
        double dt = 0.1;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double A = 1.0;
        double a = 2.1;
        double B = 0.5;
        double b = 0.8;
        double tc = 1.0;

        double C = -2.3;

        cube R = zeros<cube>(Ns, Ns, Nt);
        cube Rex = zeros<cube>(Ns, Ns, Nt);

        //Integrator
        Integrator integrator(nt, dt, ns, ds);
        vec s = integrator.coordinateVec();
        vec k = integrator.spatialFreqVec();
        vec t = integrator.timeVec();
        vec w = integrator.temporalFreqVec();


        //Cell
        DOG dog(A, a, B, b);
        TemporallyConstant Kt(tc);
        GanglionCell ganglion(&integrator, &dog, &Kt);

        //Stimulus
        double wd = w(0);
        double kx = k(1);
        double ky = k(1);
        Grating S(&integrator, {kx, ky}, wd, C);
        S.computeFourierTransform();


        //Compute analytic:
        double W = dog.fourierTransform({kx, ky}) * Kt.fourierTransform(wd);
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    Rex(i,j,l) = C * W * cos(kx*s[i] + ky*s[j] - wd * t[l]);

                }
            }
        }


        //Compute numerical
        ganglion.computeResponse(&S);
        R = ganglion.response();


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(Rex(i,j,l), R(i,j,l), 1e-12);

                }
            }
        }

    }

}

