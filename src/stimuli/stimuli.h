#ifndef STIMULI_H
#define STIMULI_H

#define PI 3.14159265359


#include <math.h>
#include <libconfig.h++>
#include <armadillo>
#include <boost/math/special_functions/bessel.hpp>

#include "../math/functions.h"

using namespace std;
using namespace libconfig;
using namespace arma;

class Stimuli
{
public:
    Stimuli(const Config *cfg);
    ~Stimuli();

    virtual double real(vec2 rVec, double t) = 0;
    virtual double complex(vec2 k, double w) = 0;

    mat real() const;
    void setReal(const mat &real);

    mat complex() const;
    void setComplex(const mat &complex);

    double w() const;

protected:
    mat m_real = zeros(2,2);
    mat m_complex = zeros(2,2);

    double m_w = 0;
    vec2 m_k = {0,0};


};

#endif // STIMULI_H
