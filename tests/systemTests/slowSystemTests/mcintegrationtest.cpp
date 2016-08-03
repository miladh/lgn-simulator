#include "mcintegrationtest.h"

MCintegrationTest::MCintegrationTest(string testLabel,
                                     string sourceFilename,
                                     double preCalls,
                                     double calls)
    : m_outputFilename("test_outputs/"+testLabel)
    , m_preCalls(preCalls)
    , m_calls(calls)
{

    cout << "==================" << testLabel <<  "==================" << endl;
    cout << "Source filename: " << sourceFilename << endl;

    string sourceFile = "systemTests/slowSystemTests/test_outputs"+sourceFilename + ".cpp";
    string hashValueFilename = "test_outputs/test_file_hash";

    // Check if file with test file ids exists------------------------------------
    if(!ifstream(hashValueFilename)){
        throw std::runtime_error("Could not open file: " + hashValueFilename);
    }


    // Read hash value of the file-----------------------------------------------
    ifstream  file_ids(hashValueFilename);
    bool hashValueIsFound = false;
    string line;
    while(getline( file_ids, line ) && !hashValueIsFound){
        if(line.find(sourceFile) != string::npos){
            m_hashValue = line.substr(0, line.find(sourceFile)) ;
            hashValueIsFound = true;
        }
    }

    // Check if file with test results exists--------------------------------
    if(!ifstream(m_outputFilename)){
        m_computeMC = true;
        cout << "output file doesn't exists: " << m_outputFilename <<", running test..." << endl;
    }

    // Check if output file has changed---------------------------------------
    else{
        cout << "output file exists: " << m_outputFilename << endl;
        ifstream  outputFile(m_outputFilename);
        bool found = false;
        string line;
        while(getline( outputFile, line ) && !found){
            if(line.find(m_hashValue) != string::npos){
                found = true;
                cout << "source file has not changed" << endl;
                outputFile.close();
                m_computeMC = false;
            }
        }
        if(!m_computeMC){
            string line;
            cout << "reading data..." << endl;
            outputFile.open(m_outputFilename);
            while(getline( outputFile, line )){
                if(line != m_hashValue){
                    m_results.push_back(stod(line));
                }
            }
            cout << "running test..." << endl;
        }else{
            cout << "source file has changed, running test with MC integration..." << endl;
            m_computeMC = true;
        }

    }
}

double integrandContainer(double *k,
                          size_t dim,
                          void *params)
{
    MCintegrationTest::Param *p = (MCintegrationTest::Param*)params;
    MCintegrationTest *c = p->a;
    return c->integrand(k, dim, params);
}

double MCintegrationTest::computeIntegral(double *xl, double *xu, void *params)
{
    double res, err;
    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_monte_function G;
    G.f = &integrandContainer;
    G.params = params;

    gsl_rng_env_setup ();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (2);
    gsl_monte_vegas_integrate (&G, xl, xu, 2, m_preCalls, r, s, &res, &err);
    //    display_results ("vegas warm-up", res, err);
    do
    {
        gsl_monte_vegas_integrate (&G, xl, xu, 2, m_calls/5, r, s, &res, &err);
        //        printf ("result = % .6f sigma = % .6f "
        //                "chisq/dof = %.1f\n", res, err,
        //                gsl_monte_vegas_chisq (s));
    }while (fabs (gsl_monte_vegas_chisq (s)) < 0.5);

    display_results ("vegas final", res, err);
    gsl_monte_vegas_free (s);
    gsl_rng_free (r);

    return res;
}


void MCintegrationTest::writeOutputFile()
{
    if(m_computeMC){
        ofstream newDataFile(m_outputFilename);
        cout << "writing new test results..." << endl;

        if (newDataFile.is_open())
        {
            newDataFile << m_hashValue << "\n";

            for(double r : m_results){
                newDataFile << std::fixed <<  std::setprecision(15) << r << endl;
            }
            newDataFile.close();
        }
    }
}


void MCintegrationTest::display_results (char *title, double result, double error)
{
    printf ("%s ==================\n", title);
    printf ("result = % .10f\n", result);
    printf ("sigma  = % .10f\n", error);
    //    printf ("exact  = % .10f\n", exact);
    //    printf ("error  = % .10f = %.2g sigma\n", result - exact,
    //            fabs (result - exact) / error);
}

