#ifndef GAUSSIAN_H
#define GAUSSIAN_H

#include "spatialkernel.h"

class Gaussian : public SpatialKernel
{
public:
    Gaussian(double weight, double spread);
    ~Gaussian();

    // SpatialKernel interface
    double spatial(vec2 rVec);
    double fourierTransform(vec2 kVec);

private:
    double m_weight = 0.0;
    double m_spread = 0;
};

#endif // GAUSSIAN_H
