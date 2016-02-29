#ifndef COMBINEDRC_H
#define COMBINEDRC_H

#include "temporalkernel.h"

namespace lgnSimulator {

class CombinedRC : public TemporalKernel
{
public:
    CombinedRC(double cenLatency, double surLatency, double delay);

    // TemporalKernel interface
    virtual double temporal(double t) const;
    virtual complex<double> fourierTransform(double w) const;

private:
    double m_cenLatency = 0.0;
    double m_surLatency = 0.0;
    double m_delay = 0.0;
};

}

lgnSimulator::CombinedRC createTemporalCombinedRCKernel(const YAML::Node &cfg);
#endif // COMBINEDRC_H
