#ifndef GAUSSIAN_H
#define GAUSSIAN_H

#include "spatialkernel.h"


using namespace std;
using namespace arma;

namespace lgnSimulator {
class SpatialGaussian : public SpatialKernel
{
public:
    SpatialGaussian(double a);

    // SpatialKernel interface
    virtual double spatial(vec2 r) const;
    virtual complex<double> fourierTransform(vec2 k) const;


private:
    double m_a = 0.0;
};
}

lgnSimulator::SpatialGaussian createSpatialGaussianKernel(const YAML::Node &cfg);

#endif // GAUSSIAN_H
