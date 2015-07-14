#ifndef SPATIALKERNELS_H
#define SPATIALKERNELS_H

#include <armadillo>

using namespace std;
using namespace arma;


class SpatialKernel
{
public:
    SpatialKernel();
    ~SpatialKernel();

    virtual double coupling(vec2 k, double w) = 0;
};

#endif // SPATIALKERNELS_H
