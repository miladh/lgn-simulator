#ifndef SPATIALLYCONSTANT_H
#define SPATIALLYCONSTANT_H

#include "spatialKernels/spatialkernel.h"


namespace edog {
class SpatiallyConstant : public SpatialKernel
{
public:
    SpatiallyConstant(double constant);

    // SpatialKernel interface
    virtual double spatial(vec2 rVec);
    virtual complex<double> fourierTransform(vec2 kVec);


private:
    double m_constant= 0.0;

};

}

edog::SpatiallyConstant createSpatiallyConstantKernel(const YAML::Node *cfg);

#endif // SPATIALLYCONSTANT_H
