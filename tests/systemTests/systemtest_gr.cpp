#include <unittest++/UnitTest++.h>

#include "neurons/ganglioncell.h"
#include "neurons/relaycell.h"

#include "stimuli/grating.h"
#include "integrator.h"

#include "../tests/systemTests/kernelsettings.h"

SUITE(SYSTEM){

    TEST(ganglion_relay){
        int ns = 5;
        int nt = 4;
        double ds = 0.1;
        double dt = 0.2;


        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        cube Rg = zeros<cube>(Ns, Ns, Nt);
        cube Rg_ex = zeros<cube>(Ns, Ns, Nt);
        cube Rr = zeros<cube>(Ns, Ns, Nt);
        cube Rr_ex = zeros<cube>(Ns, Ns, Nt);

        //Integrator
        Integrator integrator(nt, dt, ns, ds);
        vec s = integrator.coordinateVec();
        vec k = integrator.spatialFreqVec();
        vec t = integrator.timeVec();
        vec w = integrator.temporalFreqVec();

        //Stimulus
        double C = -2.3;
        double wd = w(2);
        double kx = k(1);
        double ky = k(1);
        Grating S(&integrator, {kx, ky}, wd, C);
        S.computeFourierTransform();


        //Kernels
        vector<SpatialKernel*> spatialKernels = KernelSettings::spatialKernelVector();
        vector<TemporalKernel*> temporalKernels = KernelSettings::temporalKernelVector();



        for(SpatialKernel* Fs : spatialKernels){
            for(TemporalKernel* Ft : temporalKernels){

                //ganglion cell
                GanglionCell ganglion(&integrator, Fs, Ft);
                double Wg = Fs->fourierTransform({kx, ky}) * Ft->fourierTransform(wd);

                for(SpatialKernel* Ks : spatialKernels){
                    for(TemporalKernel* Kt : temporalKernels){

                        //relay cell
                        RelayCell relay(&integrator);

                        //connect
                        relay.addGanglionCell(&ganglion, Ks, Kt);


                        //Compute analytic:
                        double Krg =  Ks->fourierTransform({kx, ky})
                                * Kt->fourierTransform(wd);
                        double Wr = Krg * Wg;

                        for(int l = 0; l < Nt; l++){
                            for(int i = 0; i < Ns; i++){
                                for(int j = 0; j < Ns; j++){
                                    Rr_ex(i,j,l) = C * Wr *
                                            cos(kx*s[i] + ky*s[j] - wd * t[l]);
                                    Rg_ex(i,j,l) = C * Wg *
                                            cos(kx*s[i] + ky*s[j] - wd * t[l]);

                                }
                            }
                        }

                        //Compute numerical
                        relay.computeResponse(&S);
                        Rr = relay.response();
                        ganglion.computeResponse(&S);
                        Rg = ganglion.response();


                        // Test
                        for(int l = 0; l < Nt; l++){
                            for(int i = 0; i < Ns; i++){
                                for(int j = 0; j < Ns; j++){
                                    CHECK_CLOSE(Rr_ex(i,j,l), Rr(i,j,l), 1e-12);
                                    CHECK_CLOSE(Rg_ex(i,j,l), Rg(i,j,l), 1e-12);

                                }
                            }
                        }

                    }//End of Kt loop
                }//End of Ks loop


            }//End of Ft loop
        }//End of Fs loop

    }

}


