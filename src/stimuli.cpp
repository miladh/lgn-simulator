#include "stimuli.h"


Stimuli::Stimuli(const Config *cfg)
{
    const Setting & root = cfg->getRoot();
    m_C = root["stimuliSettings"]["C"];
    m_d = root["stimuliSettings"]["d"];
    m_w = root["stimuliSettings"]["w"];
    m_kx = root["stimuliSettings"]["kx"];
    m_ky = root["stimuliSettings"]["ky"];

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

double Stimuli::patchGratingComplex(double kx, double ky, double w)
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


double Stimuli::delta(double x, double y) {
    if(x == y){
        return 1;
    }else{
        return 0;
    }
}



double Stimuli::w() const
{
    return m_w;
}
mat Stimuli::real() const
{
    return m_real;
}

void Stimuli::setReal(const mat &real)
{
    m_real = real;
}
mat Stimuli::complex() const
{
    return m_complex;
}

void Stimuli::setComplex(const mat &complex)
{
    m_complex = complex;
}




