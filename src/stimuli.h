#ifndef STIMULI_H
#define STIMULI_H

#define PI 3.14159265359


#include <math.h>
#include <libconfig.h++>
#include <boost/math/special_functions/bessel.hpp>



using namespace std;
using namespace libconfig;
class Stimuli
{
public:
    Stimuli(const Config *cfg);
    ~Stimuli();

    double patchGrating(double rx, double ry, double t);
    double patchGratingComplex(double kx, double ky, double w);

    double delta(double x, double y);
    double heaviside(double x);
    double secondKindBesselFunction(double x);

    double w() const;

private:
    double m_C = 0.0;
    double m_d = 0.0;
    double m_w = 0.0;
    double m_kx, m_ky = 0.0;

};

#endif // STIMULI_H
