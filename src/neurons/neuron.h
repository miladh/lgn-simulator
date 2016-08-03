#ifndef NEURON_H
#define NEURON_H

#include "stimuli/stimulus.h"
#include "integrator.h"
#include "kernels/kernel.h"
#include "staticNonlinearity/staticnonlinearity.h"

using namespace arma;
using namespace std;


namespace lgnSimulator {
class Neuron
{
public:
    Neuron(Integrator* const integrator,
           const double backgroundResponse = 0,
           const string type = "None",
           StaticNonlinearity *const staticNonlinearity = nullptr);
    ~Neuron();

    struct Input {
        Neuron* const neuron;
        const Kernel &kernel;
    };

    // Compute functions
    void computeResponse(const Stimulus *const stimulus,
                         bool recomputeFourierTransfrom=false);
    void computeResponseFourierTransform(const Stimulus *const stimulus);

    //Virtual function:
    virtual void computeImpulseResponse();
    virtual void computeImpulseResponseFourierTransform() = 0;

    // Getter member functions
    const cube &response() const;
    const cube &impulseResponse() const;
    const cx_cube &responseFourierTransform() const;
    const cx_cube &impulseResponseFourierTransform() const;

    const string type() const;

    void clearResponse();
    void clearResponseFourierTransform();
    void clearImpulseResponse();

    bool isImpulseResponseFourierTransformComputed() const;

private:
    Integrator *const m_integrator;
    StaticNonlinearity *const m_staticNonlinearity = nullptr;


protected:
    bool impulseResponseFourierTransformComputed  = false;
    bool responseFourierTransformComputed = false;
    const double m_backgroundResponse;
    const string m_type;

    cube m_response;
    cube m_impulseResponse;
    cx_cube m_responseFT;
    cx_cube m_impulseResponseFT;

    vec m_spatialVec;
    vec m_spatialFreqs;

    vec m_timeVec;
    vec m_temporalFreqs;

};
}
#endif // NEURON_H
