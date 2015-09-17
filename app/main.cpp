#include <iostream>

#include <trapezoidal.h>
#include <outputmanager.h>
#include <unistd.h>

#include "stimuli/patchgrating.h"
#include "stimuli/dogstim.h"

#include "neurons/relaycell.h"
#include "neurons/ganglioncell.h"
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

    int nSteps = root["dynamicSettings"]["nSteps"];
    double dt = root["dynamicSettings"]["dt"];

    double dogA = root["dogSettings"]["A"];
    double doga = root["dogSettings"]["a"];
    double dogB = root["dogSettings"]["B"];
    double dogb = root["dogSettings"]["b"];
    double tau_rg = root["temporalSettings"]["tau_rg"];
    double tau_rc = root["temporalSettings"]["tau_rc"];
    double delay = root["temporalSettings"]["delay"];
    double weight = 1.0;
    double spread = 0.1;


    //----------------------------------------------------------------------------

    PatchGrating S(&cfg);
    //    DOGstim S(&cfg);
    OutputManager io(&cfg);

    //Spatial kernels:
    DOG dog(dogA, doga, dogB, dogb);
    Gaussian gauss(weight, spread);
    EllipticGaussian ellipticGauss(weight, PI/4*3, 1.4, 0.1);

    //Temporal kernels:
    DecayingExponential Ktg(tau_rg, 0);
    DecayingExponential Ktc(tau_rc, delay);
    DiracDelta delta(0.0);
    DampedOscillator damped(0.425, 0.38);

    //Neurons:
    GanglionCell ganglion(&cfg, &S, &dog, &damped);
    RelayCell relay(&cfg, &S);
    CorticalCell cortical(&cfg, &S);

    vector<Neuron *> neurons;
    neurons.push_back(&ganglion);
    neurons.push_back(&relay);
    neurons.push_back(&cortical);



    relay.addGanglionCell(&ganglion,&dog, &Ktg);
    relay.addCorticalNeuron(&cortical, &ellipticGauss, &Ktg);
    cortical.addRelayCell(&relay, &dog, &delta);


    S.computeSpatiotemporal();
    S.computeFourierTransform();

    ganglion.computeResponse();
    ganglion.computeImpulseResponse();

    relay.computeResponse();
    relay.computeImpulseResponse();

//    cortical.computeResponse();
//    cortical.computeImpulseResponse();

    io.writeResponse(neurons, S);

    //    cv::Mat image;
    //    image = cv::imread("../../dog.jpg", CV_LOAD_IMAGE_GRAYSCALE);

    //    cout << double(image.at<uchar>(0,0)) << endl;
    //    cout << *reinterpret_cast<uint*>(image.data) << endl;
    //    mat arma_mat(reinterpret_cast<double*>(image.data), image.rows, image.cols );
    //    cout << arma_mat(0,0) << endl;

    return 0;
}
