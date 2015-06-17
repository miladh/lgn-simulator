#ifndef STIMULI_H
#define STIMULI_H

#define PI 3.14159265359


#include <math.h>
#include <boost/math/special_functions/bessel.hpp>


using namespace std;
class Stimuli
{
public:
    Stimuli(double C = 10.0, double d =2.5, double w =0.0,
            double kx = 10., double ky=1.);
    ~Stimuli();

    double patchGrating(double rx, double ry, double t);
    double patchGratingFT(double kx, double ky, double w);

    double delta(double x, double y);
    double heaviside(double x);
    double secondKindBesselFunction(double x);

private:
    double m_C;
    double m_d;
    double m_w;
    double m_kx, m_ky;

};

#endif // STIMULI_H
