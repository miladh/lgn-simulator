/**********************************************************************
 *  Test: response of ganglion cell, relay cell, interneuron and
 *  cortical cell with patchgrating grating stimulus:
 *
 *               Rg(r,t) = Wg(r,t) * S(r,t)
 *               Ri(r,t) = [Wg(r,t)Kig(r,t) + Kic(r,t)Kcr(r,t)*Wr(r,t)] * S(r,t)
 *               Rc(r,t) = [Wr(r,t)Kcr(r,t)] * S(r,t)
 *               Rr(r,t) = [[W(r,t)Krg(r,t) + Kri(r,t) Wg(r,t)Kig(r,t)]
 *                         /[1 - Kic(r,t)Kcr(r,t)] - Krc(r,t)Kcr(r,t)] ] * S(r,t)
 *
 *  Analytic source: monte carlo integration
 *
 * ********************************************************************/

#include <lgnSimulator.h>
#include <catch.hpp>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>


using namespace lgnSimulator;

double exact, peak;

struct Param{
    const CircleMaskGrating& S;
    const Kernel & Wg;
    const Kernel & Kig;
    const Kernel & Kic;
    const Kernel & Krg;
    const Kernel & Kri;
    const Kernel & Krc;
    const Kernel & Kcr;
};


void display_results (char *title, double result, double error)
{
    printf ("%s ==================\n", title);
    printf ("result = % .10f\n", result);
    printf ("sigma  = % .10f\n", error);
    printf ("exact  = % .10f\n", exact);
    printf ("error  = % .10f = %.2g sigma\n", result - exact,
            fabs (result - exact) / error);
}


double Wr (double *k, size_t dim, void *params)
{
    (void)(dim);
    double wd = 0.0;
    double kx = k[0];
    double ky = k[1];
    Param r = *(Param *) params;

    cx_double Wr = r.Wg.fourierTransform({kx, ky}, wd)
            *(r.Krg.fourierTransform({kx, ky}, wd)
              + r.Kri.fourierTransform({kx, ky}, wd)
              * r.Kig.fourierTransform({kx, ky}, wd))
            / (1. -r. Kic.fourierTransform({kx, ky}, wd)
               * r.Kcr.fourierTransform({kx, ky}, wd)
               * r.Kri.fourierTransform({kx, ky}, wd)
               - r.Krc.fourierTransform({kx, ky}, wd)
               * r.Kcr.fourierTransform({kx, ky}, wd));

    cx_double res =  Wr *  r.S.fourierTransformAtFrequency({kx, ky}, wd)
            /(8*core::pi*core::pi*core::pi)*peak;
    return real(res);
}

double monteCarloIntegration(Integrator integrator,
                             const CircleMaskGrating &stim,
                             const Kernel &Wg, const Kernel &Kcr,
                             const Kernel &Krg, const Kernel &Kri,const Kernel &Krc,
                             const Kernel &Kig, const Kernel &Kic,
                             size_t preCalls, size_t calls)
{
    cout << "running MC integration...." << endl;

    struct Param param = {stim, Wg, Kig, Kic, Krg, Kri, Krc, Kcr};

    double res, err;
    double xl[2] = {integrator.spatialFreqVec().min(),
                    integrator.spatialFreqVec().min()};
    double xu[2] ={integrator.spatialFreqVec().max(),
                   integrator.spatialFreqVec().max()};

    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_monte_function G;
    G.f = &Wr;
    G.params = &param;

    gsl_rng_env_setup ();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (2);
    gsl_monte_vegas_integrate (&G, xl, xu, 2, preCalls, r, s, &res, &err);
    display_results ("vegas warm-up", res, err);
    do
    {
        gsl_monte_vegas_integrate (&G, xl, xu, 2, calls/5, r, s, &res, &err);
        printf ("result = % .6f sigma = % .6f "
                "chisq/dof = %.1f\n", res, err,
                gsl_monte_vegas_chisq (s));
    }while (fabs (gsl_monte_vegas_chisq (s)) < 0.5);

    //    display_results ("vegas final", res, err);
    gsl_monte_vegas_free (s);
    gsl_rng_free (r);

    return res;

}



double runSystemTest_GRIC_pg(Integrator integrator,
                             const CircleMaskGrating &stim,
                             const Kernel &Wg, const Kernel &Kcr,
                             const Kernel &Krg, const Kernel &Kri,const Kernel &Krc,
                             const Kernel &Kig, const Kernel &Kic)
{
    //ganglion cell
    GanglionCell ganglion(&integrator, Wg);

    //relayCell cell
    RelayCell relay(&integrator);

    //cortical cell
    Interneuron interneuron(&integrator);

    //cortical cell
    CorticalCell cortical(&integrator);

    //Connect
    relay.addGanglionCell(&ganglion, Krg);
    relay.addCorticalCell(&cortical, Krc);
    relay.addInterNeuron(&interneuron, Kri);
    interneuron.addGanglionCell(&ganglion, Kig);
    interneuron.addCorticalCell(&cortical, Kic);
    cortical.addRelayCell(&relay, Kcr);

    //Compute
    ganglion.computeResponse(&stim);
    relay.computeResponse(&stim);
    interneuron.computeResponse(&stim);
    cortical.computeResponse(&stim);

    exact = relay.response()(integrator.nPointsSpatial()/2,integrator.nPointsSpatial()/2,0);

    INFO( "kd=" <<  stim.spatialFreq()
          << "  d=" << stim.maskSize()
          << "   Ns=" << integrator.nPointsSpatial());

    return exact;

}



TEST_CASE("runSystemTest_GRIC_pg_1 [slow]"){

    string sourceFilename = "systemTests/slowSystemTests/test_system_gric_patchgrating.cpp";
    string dataFilename = "test_data_GRIC_pg_1";
    string fileHash;
    vector<double> results;
    bool runTest = true;


    // Check if file with test file ids exists------------------------------------
    if(ifstream("test_file_ids")){
        cout << "Id file exists..." << endl;
        ifstream  file_ids("test_file_ids");

        bool found = false;
        string line;
        while(getline( file_ids, line ) && !found){
            if(line.find(sourceFilename) != string::npos){
                fileHash = line.substr(0, line.find(sourceFilename)) ;
                found = true;
            }
        }

        cout << fileHash << endl;
    }

    // Check if file with test results exists exists--------------------------------
    if(!ifstream(dataFilename)){
        cout << "Data file doesn't exists, creating..." << endl;
        ofstream data_file(dataFilename);
        if (data_file.is_open())
        {
            data_file << fileHash << "\n";
            data_file.close();
        }
    }else{
        ifstream  data_file(dataFilename);
        bool found = false;
        string line;
        while(getline( data_file, line ) && !found){
            if(line.find(fileHash) != string::npos){
                found = true;
                cout << "file hasn't changed" << endl;
                data_file.close();
                runTest = false;
            }
        }

        if(found){
            string line;
            vector<double> tmp;
            cout << "reading data..." << endl;
            data_file.open(dataFilename);
            while(getline( data_file, line )){
                if(line != fileHash){
                    results.push_back(stod(line));
                }
            }
        }else{
            cout << "file has changed, running test..." << endl;

        }
    }

    //-------------------------------------------------------------
    size_t preCalls = 1e4;
    size_t calls = 5e6;
    int ns = 9;
    int nt = 1;
    double dt = 1;
    double ds = 0.1;

    double w_w = 1.0;
    double w_ig = 1.0;
    double w_ic = 0.8;
    double w_rg = 1.0;
    double w_ri = -0.3;
    double w_rc = 0.8;
    double w_cr = 1.0;

    double C = 1.0;
    double orientation = 0.;
    int kId = 0.;
    int wId = 0;
    double phase = 0.0;
    vec maskSize = linspace(0.01, 13, 1);


    //integrator
    Integrator integrator(nt, dt, ns, ds);
    vec k = integrator.spatialFreqVec();
    vec w = integrator.temporalFreqVec();
    peak = integrator.temporalFreqResolution();

    DOG Ws(0.62, 1.26, 0.85);
    TemporalDelta Wt(0, dt);
    SeparableKernel W(w_w, &Ws, &Wt);


    SpatialGaussian Kig_s(1.0);
    TemporalDelta Kig_t(0, dt);
    SeparableKernel Kig(w_ig, &Kig_s, &Kig_t);

    SpatialGaussian Kic_s(1.0);
    TemporalDelta Kic_t(0, dt);
    SeparableKernel Kic(w_ic, &Kic_s, &Kic_t);


    SpatialGaussian Krg_s(0.1);
    TemporalDelta Krg_t(0, dt);
    SeparableKernel Krg(w_rg, &Krg_s, &Krg_t);

    SpatialGaussian Kri_s(0.5);
    TemporalDelta Kri_t(0, dt);
    SeparableKernel Kri(w_ri, &Kri_s, &Kri_t);

    SpatialGaussian Krc_s(0.1);
    TemporalDelta Krc_t(0, dt);
    SeparableKernel Krc(w_rc, &Krc_s, &Krc_t);

    SpatialDelta Kcr_s(ds,{0,0});
    TemporalDelta Kcr_t(0, dt);
    SeparableKernel Kcr(w_cr, &Kcr_s, &Kcr_t);

    //stimulus
    double wd = w(wId);
    double spatialFreq = k(kId);

    for(int i=0; i < maskSize.n_elem; i++){
        CircleMaskGrating stim(&integrator, spatialFreq, wd, C, phase, orientation, maskSize[i]);
        stim.computeFourierTransform();

        if(runTest){
            results.push_back(monteCarloIntegration(integrator, stim,
                                               W, Kcr,
                                               Krg, Kri, Krc,
                                               Kig, Kic,
                                               preCalls, calls));
        }
        double ftIntegrator = runSystemTest_GRIC_pg(integrator, stim,
                                                    W, Kcr,
                                                    Krg, Kri, Krc,
                                                    Kig, Kic);

        CHECK(ftIntegrator == Approx(results[i]).epsilon(1e-4));
    }

    //----------------------------------------------------------------------

    ofstream newDataFile(dataFilename);
    cout << "writing new test results..." << endl;

    if (newDataFile.is_open())
    {
        newDataFile << fileHash << "\n";

        for(double r : results){
            newDataFile << std::fixed <<  std::setprecision(15) << r;
        }
        newDataFile.close();
    }

}
