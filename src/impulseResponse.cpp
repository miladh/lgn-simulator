#include "impulseResponse.h"


ImpulseResponse::ImpulseResponse(const Config *cfg)
{
    const Setting & root = cfg->getRoot();
    m_dogA = root["dogSettings"]["A"];
    m_doga = root["dogSettings"]["a"];
    m_dogB = root["dogSettings"]["B"];
    m_dogb = root["dogSettings"]["b"];

    m_loopKernelC = root["loopKernelSettings"]["C"];
    m_loopKernelc = root["loopKernelSettings"]["c"];

    m_feedbackDelay = root["temporalSettings"]["feedbackDelay"];
    m_tau_rc = root["temporalSettings"]["tau_rc"];
    m_tau_rg = root["temporalSettings"]["tau_rg"];

}


ImpulseResponse::~ImpulseResponse()
{

}

double ImpulseResponse::edogComplex(double kx, double ky, double w)
{

    double ff = differenceOfGaussianComplex(kx,ky) * feedforwardTemporalComplex(w);
    double fb = 1.0 - loopKernel(kx,ky)* feedbackTemporalComplex(w);

    return ff/fb;
}


double ImpulseResponse::loopKernel(double kx, double ky)
{

    double k = sqrt(kx*kx + ky*ky);
    double f = m_loopKernelC * exp(-k*k * m_loopKernelc*m_loopKernelc * 0.25);

    return f;
}
mat ImpulseResponse::real() const
{
    return m_real;
}

void ImpulseResponse::setReal(const mat &real)
{
    m_real = real;

}
mat ImpulseResponse::complex() const
{
    return m_complex;
}

void ImpulseResponse::setComplex(const mat &complex)
{
    m_complex = complex;
}



double ImpulseResponse::feedforwardTemporalComplex(double w)
{
    return 1./ (1 + w*w * m_tau_rg*m_tau_rg);
}

double ImpulseResponse::feedbackTemporalComplex(double w)
{
    return cos(w* m_feedbackDelay) / (1 + w*w * m_tau_rc*m_tau_rc);
}




double ImpulseResponse::differenceOfGaussianComplex(double kx, double ky)
{


    double k = sqrt(kx*kx + ky*ky);
    double center   = m_dogA * exp(-k*k * m_doga*m_doga / 4.);
    double surround = m_dogB * exp(-k*k * m_dogb*m_dogb / 4.);

    return center - surround;
}



double ImpulseResponse::differenceOfGaussian(double rx, double ry)
{

    double r = sqrt(rx*rx + ry*ry);
    double center   = m_dogA / (m_doga*m_doga) / PI * exp(-r*r / (m_doga*m_doga));
    double surround = m_dogB / (m_dogb*m_dogb) / PI * exp(-r*r / (m_dogb*m_dogb));


    return center - surround;
}



