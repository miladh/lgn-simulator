#include <stdio.h>
#include <iostream>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <yaml-cpp/yaml.h>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

#include <lgnSimulator.h>

using namespace std;
using namespace lgnSimulator;

double exact, peak;


struct Integrand{
    Grating* S;
    DOG* K;
};

double g (double *k, size_t dim, void *params)
{
    (void)(dim);
    Integrand integrand = *(Integrand *) params;

    DOG* Wg = integrand.K;
    Grating* S = integrand.S;

    cx_double res =  Wg->fourierTransform({k[0], k[1]})
            * S->fourierTransformAtFrequency({k[0], k[1]}, 0)
            /(8*core::pi*core::pi*core::pi)*peak;

    return real(res);
}

void display_results (char *title, double result, double error)
{
    printf ("%s ==================\n", title);
    printf ("result = % .10f\n", result);
    printf ("sigma  = % .10f\n", error);
    printf ("exact  = % .10f\n", exact);
    printf ("error  = % .10f = %.2g sigma\n", result - exact,
            fabs (result - exact) / error);
}


int main(int argc, char* argv[])
{

    cout << "=====LGN Simulator Model: Spatial summation =====" << endl;

    if(argc < 2) {
        cerr << "Too few arguments." << endl;
        return 1;
    }


    //read config file-------------------------------------------------------
    YAML::Node cfg = YAML::LoadFile(argv[1]);

    //Integrator--------------------------------------------------------------
    Integrator integrator = createIntegrator(cfg["grid"]);
    peak = integrator.temporalFreqResolution();

    //Stim---------------------------------------------------------------------
    unique_ptr<Grating> S = createGratingStimulus(&integrator, cfg["stimulus"]);

    //Ganglion cell:-----------------------------------------------------------
    DOG Wg_s = createSpatialDOGKernel(cfg["ganglion"]["Wg"]);
    TemporalDelta Wg_t = createTemporalDeltaKernel(cfg["ganglion"]["Wt"]);

    SeparableKernel Wg(cfg["ganglion"]["w"].as<double>(), &Wg_s, &Wg_t);
    GanglionCell ganglion(&integrator, Wg, cfg["ganglion"]["R0"].as<double>());


    //Compute:-----------------------------------------------------------------
    S->computeFourierTransform();
    ganglion.computeImpulseResponseFourierTransform();
    ganglion.computeResponse(S.get());
    exact = ganglion.response()(integrator.nPointsSpatial()/2,integrator.nPointsSpatial()/2,0);

    //-----------------------------------------------------------------------
    double res, err;

    double xl[2] = {integrator.spatialFreqVec().min(),
                    integrator.spatialFreqVec().min()};
    double xu[2] ={integrator.spatialFreqVec().max(),
                   integrator.spatialFreqVec().max()};



    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_monte_function G;

    Integrand integrand;
    integrand.S = S.get();
    integrand.K = &Wg_s;
    G.f = &g;
    G.params = &integrand;

    size_t calls = 5e6;

    gsl_rng_env_setup ();

    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    {
        gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (2);

        gsl_monte_vegas_integrate (&G, xl, xu, 2, 3e5, r, s, &res, &err);
        display_results ("vegas warm-up", res, err);

        printf ("converging...\n");

        do
        {
            gsl_monte_vegas_integrate (&G, xl, xu, 2, calls/5, r, s, &res, &err);
            printf ("result = % .6f sigma = % .6f "
                    "chisq/dof = %.1f\n", res, err,
                    gsl_monte_vegas_chisq (s));
        }
        while (fabs (gsl_monte_vegas_chisq (s)) < 0.95);

        display_results ("vegas final", res, err);

        gsl_monte_vegas_free (s);
    }

    gsl_rng_free (r);

    return 0;
}
