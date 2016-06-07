#ifndef NONSEPARABLEDOG_H
#define NONSEPARABLEDOG_H

#include "kernel.h"
#include "spatialKernels/spatialgaussian.h"
#include "temporalKernels/doe.h"

namespace lgnSimulator{

class NonseparableDOG : public Kernel
{
public:
    NonseparableDOG(double weight, double a, double b, double c,
                    double cenLatencyAlpha, double cenLatencyBeta,
                    double surLatencyAlpha, double surLatencyBeta,
                    double delay);

    // Kernel interface
    virtual double spatiotemporal(vec2 r, double t) const;
    virtual complex<double> fourierTransform(vec2 k, double w) const;


private:
    double m_relativeStrength = 0.0;
    SpatialGaussian *m_spatialCenter;
    SpatialGaussian *m_spatialSurround;
    DOE *m_temporalCenter;
    DOE *m_temporalSurround;
};

}

lgnSimulator::NonseparableDOG createNonseparableDOGKernel(const YAML::Node &cfg);

#endif // NONSEPARABLEDOG_H
