#ifndef NEURON_H
#define NEURON_H

#include <armadillo>
#include <libconfig.h++>

#include "../lib.h"
#include "../stimuli/stimuli.h"
#include "../temporalKernels/temporalkernel.h"
#include "../spatialKernels/spatialkernel.h"

using namespace libconfig;
using namespace arma;
using namespace std;

class Neuron
{
public:
    Neuron(const Config *cfg);
    ~Neuron();

    virtual void computeResponse() = 0;

    void addFeedforwardInput(Neuron *neuron,
                             TemporalKernel *tKernel,
                             SpatialKernel *sKernel);
    void addFeedbackInput(Neuron *neuron,
                          TemporalKernel *tKernel,
                          SpatialKernel *sKernel);

    // Getter member functions
    mat response() const;
    mat impulseResponse() const;
    mat responseComplex() const;
    mat impulseResponseComplex() const;


protected:
    Stimuli *m_stim;

    mat m_response, m_responseComplex;
    mat m_impulseResponse, m_impulseResponseComplex;

    vec m_mesh;
    vec3 m_domain;

    struct Input {
        Neuron *neuron;
        TemporalKernel *temporalKernel;
        SpatialKernel *spatialKernel;
    };

    vector<Input> m_feedforwardInputs;
    vector<Input> m_feedbackInputs;




};

#endif // NEURON_H
