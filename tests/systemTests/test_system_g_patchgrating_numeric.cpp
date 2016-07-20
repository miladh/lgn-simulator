/**********************************************************************
 *  Test: ganglion cell response for patch grating, in the case where
 *        the stimulation spot and receptive-field are concentric.
 *
 *  Analytic source: Einevoll et. al (2005)
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



////Grating---------------------------------------------
//TEST_CASE("runTest_G_pg_grating_numeric"){
//    double eps = 1e-5;
//    int ns = 9;
//    double ds = 0.1;

//    vec diameters = linspace(0.0, 10, 11);
//    vec kId =linspace(0, pow(2, ns-2), 10);


//    //integrator
//    Integrator integrator(1, 0.1, ns, ds);

//    //ganglion cell
//    DOG Ws(0.62, 1.26, 0.85);
//    TemporalDelta Wt(0, 1);
//    SeparableKernel W(1, &Ws, &Wt);
//    GanglionCell ganglion(&integrator, W);

//    ganglion.computeImpulseResponse();
//    mat irf = ganglion.impulseResponse().slice(0);



//    for(double d : diameters){
//        for(double k : kId){

//            double kd = integrator.spatialFreqVec()[k];
//            CircleMaskGrating stim(&integrator, kd, 0, 0, 1, d);

//            stim.computeFourierTransform();
//            stim.computeSpatiotemporal();

//            //fft
//            ganglion.computeResponse(&stim);
//            double Rgc = ganglion.response()(integrator.nPointsSpatial()/2,
//                                             integrator.nPointsSpatial()/2,0);

//            //numeric
//            mat s = stim.spatioTemporal().slice(0);
//            double R_num = computeIntegral(s, irf)
//                    * integrator.spatialResolution()
//                    * integrator.spatialResolution()/2.;


//            SECTION(to_string(k) + "-" +  to_string(d)) {
//                cout << "fft: " << Rgc
//                     << "- bf: "  << R_num << endl;
//                REQUIRE(R_num == Approx(Rgc).epsilon(eps));
//            }
//        }
//    }

//}




