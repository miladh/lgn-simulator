#include "combinedrc.h"

using namespace lgnSimulator;

CombinedRC::CombinedRC(double cenLatency, double surLatency, double delay)
    : m_cenLatency(cenLatency)
    , m_surLatency(surLatency)
    , m_delay(delay)
{

}


double CombinedRC::temporal(double t) const
{
    if(Special::heaviside(t - m_delay) == 0) {
        return 0;
    }else{
        double center = exp(-(t - m_delay)/m_cenLatency)/m_cenLatency/m_cenLatency;
        double surround = exp(-(t - m_delay)/m_surLatency)/m_surLatency/m_surLatency;
        return center - surround;
    }
}

complex<double> CombinedRC::fourierTransform(double w) const
{
    complex<double> expFactor = exp(core::i*w*m_delay);
    complex<double> center = expFactor/m_cenLatency/m_cenLatency
                    /(1.0/ m_cenLatency - core::i*w)
                    /(1.0/ m_cenLatency - core::i*w);
    complex<double> surround = expFactor/m_surLatency/m_surLatency
                    /(1.0/ m_surLatency - core::i*w)
                    /(1.0/ m_surLatency - core::i*w);

    return center - surround;

}

CombinedRC createTemporalCombinedRCKernel(const YAML::Node &cfg)
{

    double cenLatency = cfg["cenLatency"].as<double>();
    double surLatency = cfg["surLatency"].as<double>();
    double delay = cfg["delay"].as<double>();

    return CombinedRC(cenLatency, surLatency, delay);

}
