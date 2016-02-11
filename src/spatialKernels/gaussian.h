#ifndef GAUSSIAN_H
#define GAUSSIAN_H

#include "spatialkernel.h"


using namespace std;
using namespace arma;

namespace lgnSimulator {
class Gaussian : public SpatialKernel
{
public:
    Gaussian(double A, double a);

    // SpatialKernel interface
    virtual double spatial(vec2 r);
    virtual complex<double> fourierTransform(vec2 k);


private:
    double m_A = 0.0;
    double m_a = 0.0;
};
}

lgnSimulator::Gaussian createGaussianSpatialKernel(const YAML::Node *cfg);

#endif // GAUSSIAN_H
