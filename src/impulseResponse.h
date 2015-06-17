#ifndef IMPULSERESPONSE_H
#define IMPULSERESPONSE_H

#define PI 3.14159265359

#include<iostream>
#include <math.h>
#include <libconfig.h++>



using namespace std;
using namespace libconfig;

class ImpulseResponse {



public:
    ImpulseResponse(const Config *cfg);
    ~ImpulseResponse();

    double edogComplex(double kx, double ky, double w);
    double feedforwardTemporalFT(double w);
    double feedbackTemporalFT(double w);
    double secondKindBesselFunction(double x);

    double differenceOfGaussian(double rx, double ry);

    double differenceOfGaussianFT(double kx, double ky);

    double loopKernel(double kx, double ky);

private:
    double m_dogA = 0.0;
    double m_doga = 0.0;
    double m_dogB = 0.0;
    double m_dogb = 0.0;

    double m_loopKernelC = 0.0;
    double m_loopKernelc = 0.0;


    double m_feedbackDelay = 0.0; //[s]
    double m_tau_rc = 0.0;  //[s]
    double m_tau_rg = 0.0;  //[s]

};

#endif // IMPULSERESPONSE_H
