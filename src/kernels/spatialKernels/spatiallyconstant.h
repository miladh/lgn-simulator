#ifndef SPATIALLYCONSTANT_H
#define SPATIALLYCONSTANT_H

#include "kernels/spatialKernels/spatialkernel.h"


namespace lgnSimulator {
class SpatiallyConstant : public SpatialKernel
{
public:
    SpatiallyConstant(double constant, double spatialFreqResolution);

    // SpatialKernel interface
    virtual double spatial(vec2 r) const;
    virtual complex<double> fourierTransform(vec2 k) const;


private:
    double m_constant= 0.0;
    double m_peak  = 1.0;

};

}

lgnSimulator::SpatiallyConstant createSpatiallyConstantKernel(const YAML::Node &cfg);

#endif // SPATIALLYCONSTANT_H
