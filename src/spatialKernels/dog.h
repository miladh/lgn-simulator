#ifndef DOG_H
#define DOG_H

#include "spatialkernel.h"
#include "gaussian.h"

using namespace std;
using namespace arma;

namespace lgnSimulator {
class DOG : public SpatialKernel
{
public:
    DOG(double A, double a, double B, double b);
    ~DOG();

    // SpatialKernel interface
    double spatial(vec2 r);
    complex<double> fourierTransform(vec2 k);


private:
    Gaussian *m_centre;
    Gaussian *m_surround;


};
}

lgnSimulator::DOG createSpatialDOGKernel(const YAML::Node *cfg);

#endif // DOG_H
