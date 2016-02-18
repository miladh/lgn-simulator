#include <unittest++/UnitTest++.h>

#include "neurons/ganglioncell.h"
#include "neurons/interneuron.h"
#include "neurons/relaycell.h"
#include "neurons/corticalcell.h"

#include "stimuli/grating/fullfieldgrating.h"
#include "integrator.h"

#include "../tests/systemTests/kernelsettings.h"
#include "kernels/separablekernel.h"

SUITE(SYSTEM){

    TEST(ganglion_relay_interneuron_cortical){
        int ns = 4;
        int nt = 3;

        double dt = 0.2;
        double ds = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        cube Rg = zeros<cube>(Ns, Ns, Nt);
        cube Rg_ex = zeros<cube>(Ns, Ns, Nt);
        cube Rr = zeros<cube>(Ns, Ns, Nt);
        cube Rr_ex = zeros<cube>(Ns, Ns, Nt);
        cube Ri = zeros<cube>(Ns, Ns, Nt);
        cube Ri_ex = zeros<cube>(Ns, Ns, Nt);
        cube Rc = zeros<cube>(Ns, Ns, Nt);
        cube Rc_ex = zeros<cube>(Ns, Ns, Nt);

        //Integrator
        Integrator integrator(nt, dt, ns, ds);
        vec s = integrator.spatialVec();
        vec k = integrator.spatialFreqVec();
        vec t = integrator.timeVec();
        vec w = integrator.temporalFreqVec();

        //Stimulus
        double C = -2.3;
        double wd = w(2);
        double kx = k(1);
        double ky = k(2);
        FullFieldGrating S(&integrator, {kx, ky}, wd, C);
        S.computeFourierTransform();


        //Kernels
        vector<SpatialKernel*> spatialKernels = KernelSettings::spatialKernelVector();
        vector<TemporalKernel*> temporalKernels = KernelSettings::temporalKernelVector();


        //removing some kernels to reduce computation time
        temporalKernels.erase (temporalKernels.begin()+1, temporalKernels.begin()+3);


        //--------Loop recpField of ganglion--------
        for(SpatialKernel* Fs : spatialKernels){
            for(TemporalKernel* Ft : temporalKernels){

                //ganglion cell
                SeparableKernel F(Fs, Ft);
                GanglionCell ganglion(&integrator,&F);
                complex<double> Wg = Fs->fourierTransform({kx, ky}) * Ft->fourierTransform(wd);


                //Compute ganglion response numerical
                ganglion.computeResponse(&S);
                Rg = ganglion.GanglionCell::response();

                //Compute ganglion response analytic
                for(int l = 0; l < Nt; l++){
                    for(int i = 0; i < Ns; i++){
                        for(int j = 0; j < Ns; j++){
                            Rg_ex(i,j,l) = C * abs(Wg) *
                                    cos(kx*s[i] + ky*s[j] - wd * t[l] + arg(Wg));
                        }
                    }
                }

                //--------Loop G-R connection-----------
                for(SpatialKernel* Ks_rg : spatialKernels){
                    for(TemporalKernel* Kt_rg : temporalKernels){
                        SeparableKernel K_rg(Ks_rg, Kt_rg);
                        complex<double> Krg =  Ks_rg->fourierTransform({kx, ky})
                                * Kt_rg->fourierTransform(wd);

                        //--------Loop G-I connection-----------
                        for(SpatialKernel* Ks_ig : spatialKernels){
                            for(TemporalKernel* Kt_ig : temporalKernels){
                                SeparableKernel K_ig(Ks_ig, Kt_ig);

                                complex<double> Kig =  Ks_ig->fourierTransform({kx, ky})
                                        * Kt_ig->fourierTransform(wd);

                                //--------Loop I-R connection-----------
                                for(SpatialKernel* Ks_ri : spatialKernels){
                                    for(TemporalKernel* Kt_ri : temporalKernels){
                                        SeparableKernel K_ri(Ks_ri, Kt_ri);

                                        complex<double> Kri =  Ks_ri->fourierTransform({kx, ky})
                                                * Kt_ri->fourierTransform(wd);

                                        //--------Loop R-C and C-R connection---------
                                        for(SpatialKernel* Ks_rc : spatialKernels){
                                            for(TemporalKernel* Kt_rc : temporalKernels){
                                                SeparableKernel K_rc(Ks_rc, Kt_rc);

                                                complex<double> Krc =  Ks_rc->fourierTransform({kx, ky})
                                                        * Kt_rc->fourierTransform(wd);
                                                complex<double> Kcr =  Krc;


                                                //--------Loop C-I connection-----------
                                                for(SpatialKernel* Ks_ic : spatialKernels){
                                                    for(TemporalKernel* Kt_ic : temporalKernels){
                                                        SeparableKernel K_ic(Ks_ic, Kt_ic);

                                                        complex<double> Kic =  Ks_ic->fourierTransform({kx, ky})
                                                                * Kt_ic->fourierTransform(wd);

                                                        //cells
                                                        RelayCell relay(&integrator);
                                                        Interneuron interneuron(&integrator);
                                                        CorticalCell cortical(&integrator);


                                                        //connect cells
                                                        relay.addGanglionCell(&ganglion, &K_rg);
                                                        relay.addInterNeuron(&interneuron,  &K_ri);
                                                        relay.addCorticalNeuron(&cortical,  &K_rc);

                                                        interneuron.addGanglionCell(&ganglion, &K_ig);
                                                        interneuron.addCorticalNeuron(&cortical,  &K_ic);

                                                        cortical.addRelayCell(&relay, &K_rc);


                                                        //Compute analytic:
                                                        complex<double> Wr = (Krg * Wg + Kri * Kig * Wg)
                                                                /(1. - Kri * Kic * Kcr -  Krc * Kcr);
                                                        complex<double> Wc = Kcr * Wr;
                                                        complex<double> Wi = Kig * Wg + Kic * Wc;


                                                        for(int l = 0; l < Nt; l++){
                                                            for(int i = 0; i < Ns; i++){
                                                                for(int j = 0; j < Ns; j++){
                                                                    Rr_ex(i,j,l) = C * abs(Wr) *
                                                                            cos(kx*s[i] + ky*s[j] - wd * t[l]+ arg(Wr));
                                                                    Ri_ex(i,j,l) = C * abs(Wi) *
                                                                            cos(kx*s[i] + ky*s[j] - wd * t[l] + arg(Wi));
                                                                    Rc_ex(i,j,l) = C * abs(Wc) *
                                                                            cos(kx*s[i] + ky*s[j] - wd * t[l] + arg(Wc));

                                                                }
                                                            }
                                                        }

                                                        //Compute numerical
                                                        relay.computeResponse(&S);
                                                        cortical.computeResponse(&S);
                                                        interneuron.computeResponse(&S);
                                                        Rr = relay.response();
                                                        Ri = interneuron.response();
                                                        Rc = cortical.response();



                                                        // Test
                                                        for(int l = 0; l < Nt; l++){
                                                            for(int i = 0; i < Ns; i++){
                                                                for(int j = 0; j < Ns; j++){
                                                                    CHECK_CLOSE(Rr_ex(i,j,l), Rr(i,j,l), 1e-12);
                                                                    CHECK_CLOSE(Ri_ex(i,j,l), Ri(i,j,l), 1e-12);
                                                                    CHECK_CLOSE(Rc_ex(i,j,l), Rc(i,j,l), 1e-12);

                                                                }
                                                            }
                                                        }

                                                    }//End of Kt_ic loop
                                                }//End of Ks_ic loop

                                            }//End of Kt_rc loop
                                        }//End of Ks_rc loop

                                    }//End of Kt_ri loop
                                }//End of Ks_ri loop

                            }//End of Kt_ig loop
                        }//End of Ks_ig loop

                    }//End of Kt_rg loop
                }//End of Ks_rg loop


                // Test
                for(int l = 0; l < Nt; l++){
                    for(int i = 0; i < Ns; i++){
                        for(int j = 0; j < Ns; j++){
                            CHECK_CLOSE(Rg_ex(i,j,l), Rg(i,j,l), 1e-12);

                        }
                    }
                }

            }//End of Ft loop
        }//End of Fs loop

    }

}



