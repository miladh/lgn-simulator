#include "originaledog.h"

OriginalEDOG::OriginalEDOG(const Config *cfg, Ganglion *ganglion, Stimuli *stim)
    : Relay(cfg, ganglion, stim)
{

    const Setting & root = cfg->getRoot();
    m_loopKernelC = root["loopKernelSettings"]["C"];
    m_loopKernelc = root["loopKernelSettings"]["c"];

    m_feedbackDelay = root["temporalSettings"]["feedbackDelay"];
    m_tau_rc = root["temporalSettings"]["tau_rc"];
    m_tau_rg = root["temporalSettings"]["tau_rg"];

}

OriginalEDOG::~OriginalEDOG()
{

}

double OriginalEDOG::transferFunctionComplex(vec2 kVec, double w)
{
    double k = sqrt(dot(kVec,kVec));

    double ffTemporal = 1./ (1 + w*w * m_tau_rg*m_tau_rg);

    double fbTemporal = cos(w* m_feedbackDelay) / (1 + w*w * m_tau_rc*m_tau_rc);
    double fbSpatial = m_loopKernelC * exp(-k*k * m_loopKernelc*m_loopKernelc * 0.25);

    return ffTemporal / (1 - fbTemporal*fbSpatial);
}

