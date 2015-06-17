#ifndef IMPULSERESPONSE_H
#define IMPULSERESPONSE_H

#define PI 3.14159265359

#include<iostream>
#include <math.h>



using namespace std;

class ImpulseResponse {


public:
    ImpulseResponse();
    ~ImpulseResponse();

    double edogImpulseResponseFunctionFT(double kx, double ky, double w);
    double feedforwardTemporalFT(double w);
    double feedbackTemporalFT(double w);
    double secondKindBesselFunction(double x);

    double differenceOfGaussian(double rx, double ry,
                                double A = 1., double a = 0.25,
                                double B = 0.85, double b = 0.83);

    double differenceOfGaussianFT(double kx, double ky,
                                  double A = 1., double a = 0.25,
                                  double B = 0.85, double b = 0.83);

    double loopKernel(double kx, double ky,
                      double C = 0.5, double c = 0.83);

private:
    double m_feedbackDelay = 0.1; //[s]
    double m_tau_rc = 0.42;  //[s]
    double m_tau_rg = 0.26;  //[s]

};

#endif // IMPULSERESPONSE_H
