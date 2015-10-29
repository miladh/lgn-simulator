#include <unittest++/UnitTest++.h>

#include "neurons/ganglioncell.h"
#include "neurons/relaycell.h"

#include "stimuli/grating.h"

#include "integrator.h"

#include "spatialKernels/dog.h"
#include "spatialKernels/gaussian.h"
#include "spatialKernels/ellipticgaussian.h"

#include "temporalKernels/decayingexponential.h"
#include "temporalKernels/temporallyconstant.h"


SUITE(SYSTEM){

    TEST(ganglion_relay){
        int ns = 5;
        int nt = 5;
        double ds = 0.1;
        double dt = 0.1;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        //Stim
        double C = -2.3;

        //Decaying exp temporal kernel
        double tau = 1.0;
        double delay = 0.0;

        //DOG spatial kernel
        double A = 1.0;
        double a = 2.1;
        double B = 0.5;
        double b = 0.8;

        //Gauss spatial kernel
        double weight = 0.5;
        double spread = 1.1;

        //constant temporal kernel:
        double constant = 3.1;

        cube R = zeros<cube>(Ns, Ns, Nt);
        cube Rex = zeros<cube>(Ns, Ns, Nt);

        //Integrator
        Integrator integrator(nt, dt, ns, ds);
        vec s = integrator.coordinateVec();
        vec k = integrator.spatialFreqVec();
        vec t = integrator.timeVec();
        vec w = integrator.temporalFreqVec();


        //Ganglion cell
        DOG dog(A, a, B, b);
        DecayingExponential Kt(tau, delay);
        GanglionCell ganglion(&integrator, &dog, &Kt);

        //Relay cell
        RelayCell relay(&integrator);

        //connect
        Gaussian gauss(weight, spread);
        TemporallyConstant tempConst(constant);
        relay.addGanglionCell(&ganglion, &gauss, &tempConst);

        //Stimulus
        double wd = w(0);
        double kx = k(1);
        double ky = k(1);
        Grating S(&integrator, {kx, ky}, wd, C);
        S.computeFourierTransform();


        //Compute analytic:
        double Wg =   dog.fourierTransform({kx, ky}) * Kt.fourierTransform(wd);
        double Krg =  gauss.fourierTransform({kx, ky})* tempConst.fourierTransform(wd);
        double Wr = Krg * Wg;

        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    Rex(i,j,l) = C * Wr * cos(kx*s[i] + ky*s[j] - wd * t[l]);

                }
            }
        }


        //Compute numerical
        relay.computeResponse(&S);
        R = relay.response();

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


