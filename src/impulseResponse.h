#ifndef IMPULSERESPONSE_H
#define IMPULSERESPONSE_H

#define PI 3.14159265359

#include<iostream>
#include <math.h>
#include <libconfig.h++>
#include <armadillo>


using namespace std;
using namespace libconfig;
using namespace arma;

class ImpulseResponse {



public:
    ImpulseResponse(const Config *cfg);
    ~ImpulseResponse();

    double edogComplex(double kx, double ky, double w);
    double feedforwardTemporalComplex(double w);
    double feedbackTemporalComplex(double w);
    double secondKindBesselFunction(double x);

    double differenceOfGaussian(double rx, double ry);

    double differenceOfGaussianComplex(double kx, double ky);

    double loopKernel(double kx, double ky);

    mat real() const;
    void setReal(const mat &real);

    mat complex() const;
    void setComplex(const mat &complex);

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



    mat m_real = zeros(2,2);
    mat m_complex = zeros(2,2);





};

#endif // IMPULSERESPONSE_H
