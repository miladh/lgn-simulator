#include "impulseResponse.h"


ImpulseResponse::ImpulseResponse()
{

}

ImpulseResponse::~ImpulseResponse()
{

}

double ImpulseResponse::edogImpulseResponseFunction(double kx, double ky, double w)
{

    double ff = differenceOfGaussianFT(kx,ky) * feedforwardTemporalFT(w);
    double fb = 1.0 - loopKernel(kx,ky)* feedbackTemporalFT(w);

    return ff/fb;
}


double ImpulseResponse::loopKernel(double kx, double ky, double C, double c)
{

    double k = sqrt(kx*kx + ky*ky);
    double f = C * exp(-k*k * c*c * 0.25);

    return f;
}

double ImpulseResponse::feedforwardTemporalFT(double w)
{
    return 1.0;
}

double ImpulseResponse::feedbackTemporalFT(double w)
{
    return 1.0;
}




double ImpulseResponse::differenceOfGaussianFT(double kx, double ky,
                                               double A, double a,
                                               double B, double b)
{


    double k = sqrt(kx*kx + ky*ky);
    double center   = A * exp(-k*k * a*a / 4.);
    double surround = B * exp(-k*k * b*b / 4.);

    return center - surround;
}



double ImpulseResponse::differenceOfGaussian(double rx, double ry,
                                             double A, double a,
                                             double B, double b)
{

    double r = sqrt(rx*rx + ry*ry);
    double center   = A / (a*a) / PI * exp(-r*r / (a*a));
    double surround = B / (b*b) / PI * exp(-r*r / (b*b));


    return center - surround;
}



