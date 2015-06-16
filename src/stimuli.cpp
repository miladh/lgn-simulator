#include "stimuli.h"

Stimuli::Stimuli()
{

}

Stimuli::~Stimuli()
{

}


double Stimuli::patchGrating(double rx, double ry, double t,
                             double kx, double ky, double w,
                             double C, double d)
{
    double r = sqrt(rx*rx + ry*ry);
    double s = C * (1 -  heaviside(r - d * 0.5)) * cos(kx*rx + ky*ry - w*t);

    return s;
}

double Stimuli::patchGratingFT(double kx, double ky,double w,
                               double kpgx, double kpgy, double wpg,
                               double C, double d)
{
    double arg = sqrt(kx - kpgx + ky - kpgy);
    double s = secondKindBesselFunction(arg * d * 0.5) * delta(w, wpg);

    return s*2*PI*PI*d*C/arg;

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



