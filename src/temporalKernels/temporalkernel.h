#ifndef TEMPORALKERNELS_H
#define TEMPORALKERNELS_H

#include <armadillo>
#include <libconfig.h++>

#include "integrator.h"
#include "../math/functions.h"

using namespace std;
using namespace arma;
using namespace libconfig;

class TemporalKernel
{
public:
    TemporalKernel();
    ~TemporalKernel();

    virtual double temporal(double t) = 0;
    virtual double fourierTransform(double w) = 0;

};

#endif // TEMPORALKERNELS_H
