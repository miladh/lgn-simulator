#include <iostream>

#include <outputmanager.h>
#include <unistd.h>

#include "stimuli/patchgrating.h"
#include "stimuli/grating.h"

#include "integrator.h"
#include "integratorsettings.h"

#include "neurons/ganglioncell.h"
#include "neurons/interneuron.h"
#include "neurons/relaycell.h"
#include "neurons/corticalcell.h"

#include "spatialKernels/dog.h"
#include "spatialKernels/gaussian.h"
#include "spatialKernels/ellipticgaussian.h"

#include "temporalKernels/decayingexponential.h"
#include "temporalKernels/diracDelta.h"
#include "temporalKernels/dampedoscillator.h"


#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>


using namespace std;

int main()
{

    cout << "=====Extended-DOG Model=====" << endl;

    //read config file---------------------------------------------------------------
    Config cfg;
    cfg.readFile("../../eDOG/app/config.cfg");
    const Setting & root = cfg.getRoot();


    double dogA = root["dogSettings"]["A"];
    double doga = root["dogSettings"]["a"];
    double dogB = root["dogSettings"]["B"];
    double dogb = root["dogSettings"]["b"];
    double tau_rg = root["temporalSettings"]["tau_rg"];
    double tau_rc = root["temporalSettings"]["tau_rc"];
    double delay = root["temporalSettings"]["delay"];
    double weight = 0.5;
    double spread = 1.1;

    int ns = root["integratorSettings"]["ns"];
    int nt = root["integratorSettings"]["nt"];
    double dt = root["integratorSettings"]["dt"];
    double ds = root["integratorSettings"]["ds"];

    double contrast = root["stimuliSettings"]["C"];
    double spotDiameter = root["stimuliSettings"]["d"];


    //----------------------------------------------------------------------------

    IntegratorSettings integratorSettings(nt, dt, ns, ds);
    Integrator integrator(&integratorSettings);

    vec k = integrator.spatialFreqVec();
    vec w = integrator.temporalFreqVec();
    double wd = w(w.n_elem/2+2);
    double kx = k(k.n_elem/2+3);
    double ky = k(k.n_elem/2);
//    Grating S(integrator, {kx, ky}, wd, contrast);
    PatchGrating S(integrator, {kx, ky}, wd, contrast, spotDiameter);
    OutputManager io(&cfg);

    //Spatial kernels:
    DOG dog(dogA, doga, dogB, dogb);
//    Gaussian gauss(weight, spread);
//    EllipticGaussian ellipticGauss(weight, PI/4*3, 0.1, 1.0);

    //Temporal kernels:
    DecayingExponential Ktg(tau_rg, 1);
//    DecayingExponential Ktc(tau_rc, delay);
//    DiracDelta delta(0.0);
//    DampedOscillator damped(2.425, 1.38);

    //Neurons:
    GanglionCell ganglion(&cfg, &S, integrator, &dog, &Ktg);
//    RelayCell relay(&cfg, &S, integrator);
//    Interneuron interneuron(&cfg, &S, integrator);
//    CorticalCell cortical(&cfg, &S, integrator);

    vector<Neuron *> neurons;
    neurons.push_back(&ganglion);
//    neurons.push_back(&relay);
//    neurons.push_back(&interneuron);
//    neurons.push_back(&cortical);


//    interneuron.addGanglionCell(&ganglion,&dog, &Ktg);
//    interneuron.addCorticalNeuron(&cortical, &ellipticGauss, &Ktc);

//    relay.addGanglionCell(&ganglion,&gauss, &Ktg);
//    relay.addInterNeuron(&interneuron,&dog, &damped);
//    relay.addCorticalNeuron(&cortical, &ellipticGauss, &Ktc);

//    cortical.addRelayCell(&relay, &dog, &Ktg);


    //Compute:
    S.computeSpatiotemporal();
    S.computeFourierTransform();

    ganglion.computeResponse();
    ganglion.computeImpulseResponse();

//    interneuron.computeResponse();
//    interneuron.computeImpulseResponse();

//    relay.computeResponse();
//    relay.computeImpulseResponse();

//    cortical.computeResponse();
//    cortical.computeImpulseResponse();

    io.writeResponse(neurons, S);


    return 0;
}

