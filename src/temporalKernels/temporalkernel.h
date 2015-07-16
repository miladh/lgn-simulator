#ifndef TEMPORALKERNELS_H
#define TEMPORALKERNELS_H

#include <armadillo>
#include "../math/functions.h"

using namespace arma;
using namespace std;

class TemporalKernel
{
public:
    TemporalKernel();
    ~TemporalKernel();

    virtual double real(double t) = 0;
    virtual double complex(double w) = 0;

protected:
    Functions m_math;

};

#endif // TEMPORALKERNELS_H
