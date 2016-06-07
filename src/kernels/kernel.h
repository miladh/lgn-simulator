#ifndef KERNEL_H
#define KERNEL_H

#include <armadillo>
#include <yaml-cpp/yaml.h>
#include "helper/special.h"

using namespace std;
using namespace arma;

namespace lgnSimulator {
class Kernel
{
public:
    Kernel(double weight);

    virtual double spatiotemporal(vec2 r, double t) const = 0;
    virtual complex<double> fourierTransform(vec2 k, double w) const = 0;

    double weight() const;

protected:
    double m_weight = 0.0;
};

}
#endif // KERNEL_H
