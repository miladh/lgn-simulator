#ifndef TEMPORALKERNELS_H
#define TEMPORALKERNELS_H

#include <armadillo>

using namespace arma;
using namespace std;

class TemporalKernel
{
public:
    TemporalKernel();
    ~TemporalKernel();

    virtual double coupling(vec2 k, double w) = 0;
};

#endif // TEMPORALKERNELS_H
