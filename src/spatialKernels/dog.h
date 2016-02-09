#ifndef DOG_H
#define DOG_H

#include "spatialkernel.h"

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
    double m_A = 0.0;
    double m_a = 0.0;
    double m_B = 0.0;
    double m_b = 0.0;


};
}

lgnSimulator::DOG createDOGSpatialKernel(const YAML::Node *cfg);

#endif // DOG_H
