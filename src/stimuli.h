#ifndef STIMULI_H
#define STIMULI_H

#define PI 3.14159265359


#include <math.h>
#include <libconfig.h++>
#include <armadillo>
#include <boost/math/special_functions/bessel.hpp>



using namespace std;
using namespace libconfig;
using namespace arma;

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

    mat real() const;
    void setReal(const mat &real);

    mat complex() const;
    void setComplex(const mat &complex);

private:
    double m_C = 0.0;
    double m_d = 0.0;
    double m_w = 0.0;
    double m_kx, m_ky = 0.0;

    mat m_real = zeros(2,2);
    mat m_complex = zeros(2,2);


};

#endif // STIMULI_H
