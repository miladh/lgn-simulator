#ifndef SPATIALKERNELS_H
#define SPATIALKERNELS_H

#include <armadillo>
#include <yaml-cpp/yaml.h>
#include "helper/specialfunctions.h"
#include "helper/helperconstants.h"

using namespace std;
using namespace arma;

namespace lgnSimulator {
class SpatialKernel
{
public:
    SpatialKernel();
    ~SpatialKernel();

    virtual double spatial(vec2 r) const = 0;
    virtual complex<double> fourierTransform(vec2 k) const = 0;

};
}
#endif // SPATIALKERNELS_H
