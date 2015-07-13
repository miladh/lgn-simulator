#ifndef ORIGINALEDOG_H
#define ORIGINALEDOG_H


#include "relay.h"

class OriginalEDOG : public Relay
{
public:
    OriginalEDOG(const Config *cfg, Ganglion *ganglion, Stimuli * stim);
    ~OriginalEDOG();

    // Relay interface
public:
    double transferFunctionComplex(vec2 kVec, double w);


private:
    double m_loopKernelC = 0.0;
    double m_loopKernelc = 0.0;


    double m_feedbackDelay = 0.0; //[s]
    double m_tau_rc = 0.0;  //[s]
    double m_tau_rg = 0.0;  //[s]
};

#endif // ORIGINALEDOG_H
