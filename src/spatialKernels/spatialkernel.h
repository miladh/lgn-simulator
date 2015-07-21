#ifndef SPATIALKERNELS_H
#define SPATIALKERNELS_H

#include <armadillo>
#include "math/functions.h"

using namespace std;
using namespace arma;


class SpatialKernel
{
public:
    SpatialKernel();
    ~SpatialKernel();

    virtual double real(vec2 rVec) = 0;
    virtual double complex(vec2 kVec) = 0;

};

#endif // SPATIALKERNELS_H
