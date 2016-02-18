#ifndef KERNEL_H
#define KERNEL_H

#include <armadillo>
#include <yaml-cpp/yaml.h>
#include "math/functions.h"

using namespace std;
using namespace arma;

namespace lgnSimulator {
class Kernel
{
public:
    Kernel();

    virtual double spatiotemporal(vec2 r, double t) const = 0;
    virtual complex<double> fourierTransform(vec2 k, double w) const = 0;

};

}
#endif // KERNEL_H
