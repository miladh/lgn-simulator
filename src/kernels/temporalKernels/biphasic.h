#ifndef BIPHASIC_H
#define BIPHASIC_H


#include "kernels/temporalKernels/temporalkernel.h"

namespace lgnSimulator {

class Biphasic : public TemporalKernel
{
public:
    Biphasic(double phaseDuration, double dampingFactor, double delay);
    ~Biphasic();

    // TemporalKernel interface
    double temporal(double t) const;
    complex<double> fourierTransform(double w) const;


private:
    double m_phaseDuration = 1.0;
    double m_dampingFactor = 0.0;
    double m_delay = 0.0;
};

}
lgnSimulator::Biphasic createTemporalBiphasicKernel(const YAML::Node &cfg);

#endif // BIPHASIC_H
