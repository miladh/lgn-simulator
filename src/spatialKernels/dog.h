#ifndef DOG_H
#define DOG_H

#include <armadillo>
#include "spatialkernel.h"

using namespace std;
using namespace arma;


class DOG : public SpatialKernel
{
public:
    DOG(double A, double a, double B, double b);
    ~DOG();

    // SpatialKernel interface
    double spatial(vec2 r);
    double fourierTransform(vec2 k);


private:
    double m_A = 0.0;
    double m_a = 0.0;
    double m_B = 0.0;
    double m_b = 0.0;


};

DOG createDOGSpatialKernel(const YAML::Node *cfg);

#endif // DOG_H
