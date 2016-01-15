#ifndef SPATIALKERNELS_H
#define SPATIALKERNELS_H

#include <armadillo>
#include <yaml-cpp/yaml.h>

#include "math/functions.h"

using namespace std;
using namespace arma;

namespace edog {
class SpatialKernel
{
public:
    SpatialKernel();
    ~SpatialKernel();

    virtual double spatial(vec2 rVec) = 0;
    virtual double fourierTransform(vec2 kVec) = 0;

};
}
#endif // SPATIALKERNELS_H
