/**********************************************************************
 *  Test: ganglion cell response for patch grating, in the case where
 *        the stimulation spot and receptive-field are concentric.
 *
 *  Analytic source: numerical integration
 *
 * ********************************************************************/


#include <lgnSimulator.h>
#include <catch.hpp>

using namespace lgnSimulator;

double computeIntegral(mat s, mat irf){
    double res = 0;
    for(int i=0; i < int(s.n_rows)-1; i++){
        for(int j=0; j < int(s.n_cols)-1; j++){
            res+= (s(i+1,j+1) * irf(i+1,j+1) + s(i,j) * irf(i,j));
        }
    }
    return res;
}


//Grating---------------------------------------------
TEST_CASE("runTest_G_pg_grating_numeric"){
    double eps = 1e-3;
    int ns = 10;
    double ds = 0.01;

    vec diameters = linspace(1.0, 10, 10);
    vec kId =linspace(0, pow(2, ns-2), 10);


    //integrator
    Integrator integrator(1, 0.1, ns, ds);

    //kernel
    DOG Ws(0.62, 1.26, 0.85);
    TemporalDelta Wt(0, 1);
    SeparableKernel W(1, &Ws, &Wt);

    int counter=-1;
    for(double d : diameters){
        for(double k : kId){
            counter+=1;

            //ganglion
            GanglionCell ganglion(&integrator, W);
            ganglion.computeImpulseResponse();
            mat irf = ganglion.impulseResponse().slice(0);

            //stim
            double kd = integrator.spatialFreqVec()[k];
            CircleMaskGrating stim(&integrator, kd, 0, 0, 1, d);
            stim.computeFourierTransform();
            stim.computeSpatiotemporal();

            //fft
            ganglion.computeResponse(&stim);
            double Rgc = ganglion.response()(integrator.nPointsSpatial()/2,
                                             integrator.nPointsSpatial()/2,0);

            //numeric
            mat s = stim.spatioTemporal().slice(0);
            double R_num = computeIntegral(s, irf)
                    * integrator.spatialResolution()
                    * integrator.spatialResolution()/2.;

//            cout << setprecision(10)
//                 <<"Test: " + to_string(counter)
//                 << " bf: "  << R_num
//                 << " - fft: " << Rgc
//                 << endl;

                INFO( "kd=" <<  stim.spatialFreq() << "  d=" << stim.maskSize());
                CHECK(R_num == Approx(Rgc).epsilon(eps));

        }
    }

}




