#include "stimuli.h"


Stimuli::Stimuli(double C, double d, double w, double kx, double ky)
    : m_C(C)
    , m_d(d)
    , m_w(w)
    , m_kx(kx)
    , m_ky(ky)
{

}

Stimuli::~Stimuli()
{

}

double Stimuli::patchGrating(double rx, double ry, double t)
{
    double r = sqrt(rx*rx + ry*ry);
    double s = m_C*(1 - heaviside(r - m_d * 0.5)) * cos(m_kx*rx + m_ky*ry - m_w*t);

    return s;
}

double Stimuli::patchGratingFT(double kx, double ky, double w)
{
    double arg = sqrt((kx - m_kx) * (kx - m_kx)  + (ky - m_ky) * (ky - m_ky));
    double s = secondKindBesselFunction(arg * m_d * 0.5) * delta(w, m_w);

    return s*2*PI*PI*m_d*m_C/arg;

}

double Stimuli::heaviside(double x)
{
    if (x < 0){
        return 0;
    }else{
        return 1;
    }
}

double Stimuli::secondKindBesselFunction(double x)
{
    double j =  boost::math::cyl_bessel_j(1, x);

    return j;
}

int Stimuli::delta(int x, int y) {
    return (x == y);
}



