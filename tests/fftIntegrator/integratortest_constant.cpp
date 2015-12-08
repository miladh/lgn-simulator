#include <unittest++/UnitTest++.h>
#include <armadillo>
#include <iostream>

#include "integrator.h"

using namespace std;
using namespace arma;



SUITE(INTEGRATOR){


    TEST(constant_0){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -60.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_1){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -60.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_2){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -60.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_3){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -60.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_4){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -60.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_5){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -60.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_6){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -60.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_7){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -60.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_8){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -60.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_9){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -60.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_10){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -60.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_11){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -60.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_12){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -60.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_13){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -60.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_14){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -60.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_15){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -60.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_16){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -60.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_17){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -60.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_18){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -60.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_19){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -60.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_20){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -60.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_21){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -60.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_22){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -60.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_23){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -60.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_24){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -60.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_25){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -60.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_26){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -60.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_27){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -60.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_28){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -60.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_29){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -60.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_30){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -60.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_31){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -60.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_32){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -60.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_33){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -60.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_34){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -60.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_35){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -60.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_36){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -0.002;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_37){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -0.002;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_38){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -0.002;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_39){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -0.002;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_40){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -0.002;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_41){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -0.002;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_42){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -0.002;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_43){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -0.002;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_44){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -0.002;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_45){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -0.002;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_46){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -0.002;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_47){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -0.002;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_48){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -0.002;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_49){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -0.002;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_50){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -0.002;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_51){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -0.002;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_52){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -0.002;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_53){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -0.002;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_54){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -0.002;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_55){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -0.002;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_56){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -0.002;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_57){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -0.002;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_58){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -0.002;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_59){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -0.002;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_60){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -0.002;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_61){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -0.002;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_62){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -0.002;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_63){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -0.002;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_64){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -0.002;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_65){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -0.002;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_66){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -0.002;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_67){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -0.002;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_68){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -0.002;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_69){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -0.002;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_70){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -0.002;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_71){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = -0.002;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_72){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 0.0;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_73){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 0.0;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_74){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 0.0;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_75){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 0.0;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_76){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 0.0;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_77){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 0.0;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_78){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 0.0;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_79){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 0.0;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_80){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 0.0;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_81){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 0.0;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_82){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 0.0;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_83){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 0.0;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_84){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 0.0;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_85){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 0.0;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_86){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 0.0;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_87){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 0.0;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_88){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 0.0;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_89){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 0.0;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_90){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 0.0;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_91){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 0.0;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_92){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 0.0;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_93){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 0.0;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_94){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 0.0;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_95){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 0.0;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_96){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 0.0;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_97){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 0.0;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_98){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 0.0;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_99){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 0.0;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_100){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 0.0;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_101){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 0.0;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_102){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 0.0;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_103){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 0.0;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_104){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 0.0;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_105){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 0.0;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_106){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 0.0;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_107){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 0.0;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_108){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 4.6;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_109){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 4.6;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_110){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 4.6;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_111){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 4.6;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_112){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 4.6;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_113){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 4.6;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_114){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 4.6;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_115){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 4.6;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_116){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 4.6;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_117){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 4.6;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_118){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 4.6;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_119){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 4.6;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_120){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 4.6;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_121){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 4.6;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_122){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 4.6;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_123){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 4.6;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_124){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 4.6;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_125){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 4.6;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_126){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 4.6;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_127){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 4.6;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_128){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 4.6;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_129){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 4.6;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_130){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 4.6;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_131){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 4.6;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_132){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 4.6;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_133){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 4.6;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_134){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 4.6;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_135){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 4.6;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_136){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 4.6;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_137){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 4.6;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_138){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 4.6;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_139){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 4.6;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_140){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 4.6;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_141){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 4.6;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_142){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 4.6;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_143){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 4.6;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_144){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 32.1;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_145){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 32.1;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_146){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 32.1;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_147){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 32.1;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_148){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 32.1;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_149){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 32.1;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_150){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 32.1;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_151){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 32.1;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_152){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 32.1;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_153){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 32.1;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_154){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 32.1;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_155){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 32.1;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_156){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 32.1;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_157){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 32.1;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_158){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 32.1;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_159){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 32.1;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_160){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 32.1;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_161){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 32.1;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_162){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 32.1;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_163){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 32.1;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_164){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 32.1;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_165){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 32.1;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_166){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 32.1;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_167){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 32.1;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_168){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 32.1;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_169){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 32.1;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_170){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 32.1;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_171){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 32.1;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_172){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 32.1;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_173){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 32.1;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_174){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 32.1;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_175){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 32.1;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_176){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 32.1;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_177){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 32.1;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_178){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 32.1;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_179){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 32.1;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_180){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 1000.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_181){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 1000.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_182){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 1000.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_183){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 1000.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_184){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 1000.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_185){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 1000.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_186){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 1000.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_187){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 1000.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_188){
        int ns = 1;
        int nt = 1;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 1000.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_189){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 1000.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_190){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 1000.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_191){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 1000.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_192){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 1000.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_193){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 1000.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_194){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 1000.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_195){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 1000.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_196){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 1000.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_197){
        int ns = 1;
        int nt = 2;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 1000.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_198){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 1000.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_199){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 1000.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_200){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 1000.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_201){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 1000.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_202){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 1000.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_203){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 1000.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_204){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 1000.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_205){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 1000.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_206){
        int ns = 2;
        int nt = 1;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 1000.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_207){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 1000.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_208){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 1000.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_209){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 1000.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_210){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 1000.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_211){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 1000.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_212){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 1000.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_213){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.01;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 1000.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_214){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.255;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 1000.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }
    TEST(constant_215){
        int ns = 2;
        int nt = 2;
        
        double dt = 0.5;

        int Ns = pow(2,ns);
        int Nt = pow(2,nt);

        double c = 1000.5;

        Integrator integrator(nt, dt, ns);

        vec k = integrator.spatialFreqVec();
        vec w = integrator.temporalFreqVec();

        cx_cube g = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube G = zeros<cx_cube>(Ns, Ns, Nt);
        cx_cube f = zeros<cx_cube>(Ns, Ns, Nt);


        //Spatiotemporal signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    g(i,j,l) = c;
                }
            }
        }

        //fourier signal
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    f(i,j,l) = Functions::delta(k[i], 0)
                            * Functions::delta(k[j], 0)
                            * Functions::delta(w[l], 0);
                }
            }
        }

        f *= c*8*PI*PI*PI;
        f /= integrator.spatialFreqResolution()
                * integrator.spatialFreqResolution()
                * integrator.temporalFreqResolution();

        // Backward
        G = integrator.backwardFFT(f);


        // Test
        for(int l = 0; l < Nt; l++){
            for(int i = 0; i < Ns; i++){
                for(int j = 0; j < Ns; j++){
                    CHECK_CLOSE(real(g(i,j,l)),
                                real(G(i,j,l)), 1e-9);
                }
            }
        }
    }

}



