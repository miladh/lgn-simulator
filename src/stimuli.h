#ifndef STIMULI_H
#define STIMULI_H

#define PI 3.14159265359

#include <math.h>
#include <boost/math/special_functions/bessel.hpp>

class Stimuli
{
public:
    Stimuli();
    ~Stimuli();

    double patchGrating(double rx, double ry, double t,
                        double kx = 10., double ky=1., double w=0.1,
                        double C = 10.0, double d =2.5 );

    double patchGratingFT(double kx, double ky, double w,
                          double kpgx = 10., double kpgy=1., double wpg =0.1,
                          double C = 10.0, double d =2.5 );

    int delta(int x, int y);
    double heaviside(double x);
    double secondKindBesselFunction(double x);
};

#endif // STIMULI_H
