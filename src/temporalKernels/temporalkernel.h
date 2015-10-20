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

    virtual double temporal(double t) = 0;
    virtual double fourierTransform(double w) = 0;


};

#endif // TEMPORALKERNELS_H
