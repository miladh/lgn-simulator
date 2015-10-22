#ifndef SPATIALKERNELS_H
#define SPATIALKERNELS_H

#include <armadillo>
#include <libconfig.h++>

#include "math/functions.h"

using namespace std;
using namespace arma;
using namespace libconfig;


class SpatialKernel
{
public:
    SpatialKernel();
    ~SpatialKernel();

    virtual double spatial(vec2 rVec) = 0;
    virtual double fourierTransform(vec2 kVec) = 0;

};

#endif // SPATIALKERNELS_H
