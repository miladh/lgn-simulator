#include "spatiallyconstant.h"

using namespace edog;

SpatiallyConstant::SpatiallyConstant(double constant)
    : m_constant(constant)
{

}

double edog::SpatiallyConstant::spatial(vec2 rVec)
{
    (void)rVec;
    return m_constant;
}

complex<double> edog::SpatiallyConstant::fourierTransform(vec2 kVec)
{
    return m_constant * Functions::delta(kVec(0),0) * Functions::delta(kVec(1),0);
}



edog::SpatiallyConstant createSpatiallyConstantSpatialKernel(const YAML::Node *cfg)
{
    double constant = (*cfg)["constant_spatiallyDelta"].as<double>();

    return SpatiallyConstant(constant);

}
