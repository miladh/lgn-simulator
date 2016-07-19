#ifndef DOG_H
#define DOG_H

#include "spatialkernel.h"
#include "spatialgaussian.h"

using namespace std;
using namespace arma;

namespace lgnSimulator {
class DOG : public SpatialKernel
{
public:
    explicit DOG(double a, double b, double c);
    DOG(const DOG &dog);
    ~DOG();

    // SpatialKernel interface
    double spatial(vec2 r) const;
    complex<double> fourierTransform(vec2 k) const;


    double relativeStrength() const;

private:
    double m_relativeStrength = 1.0;
    SpatialGaussian* m_center = nullptr;
    SpatialGaussian* m_surround = nullptr;



};
}

lgnSimulator::DOG createSpatialDOGKernel(const YAML::Node &cfg);

#endif // DOG_H
