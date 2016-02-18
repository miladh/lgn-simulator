#ifndef SPATIALKERNELS_H
#define SPATIALKERNELS_H

#include <armadillo>
#include <yaml-cpp/yaml.h>
#include "math/functions.h"

using namespace std;
using namespace arma;

namespace lgnSimulator {
class SpatialKernel
{
public:
    SpatialKernel();
    ~SpatialKernel();

    virtual double spatial(vec2 r) = 0;
    virtual complex<double> fourierTransform(vec2 k) = 0;

};
}
#endif // SPATIALKERNELS_H
