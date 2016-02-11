#include "spatiallyconstant.h"

using namespace lgnSimulator;

SpatiallyConstant::SpatiallyConstant(double constant)
    : m_constant(constant)
{

}

double lgnSimulator::SpatiallyConstant::spatial(vec2 r)
{
    (void)r;
    return m_constant;
}

complex<double> lgnSimulator::SpatiallyConstant::fourierTransform(vec2 k)
{
    return m_constant * Functions::delta(k(0),0) * Functions::delta(k(1),0);
}



lgnSimulator::SpatiallyConstant createSpatiallyConstantKernel(const YAML::Node *cfg)
{
    double constant = (*cfg)["constant"].as<double>();

    return SpatiallyConstant(constant);

}