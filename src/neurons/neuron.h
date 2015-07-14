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

    struct Input {
        Neuron *neuron;
        TemporalKernel *temporalKernel;
        SpatialKernel *spatialKernel;
    };

    virtual void computeResponse() = 0;
    virtual double impulseResponseComplex(vec2 kVec, double w) = 0;

    void addGanglionCell(Neuron *neuron,
                          TemporalKernel *tKernel,
                          SpatialKernel *sKernel);
    void addInterNeuron(Neuron *neuron,
                             TemporalKernel *tKernel,
                             SpatialKernel *sKernel);
    void addCorticalNeuron(Neuron *neuron,
                          TemporalKernel *tKernel,
                          SpatialKernel *sKernel);

    // Getter member functions
    mat response() const;
    mat impulseResponse() const;
    mat responseComplex() const;
    mat impulseResponseComplex() const;


    vector<Input> ganglionCells() const;
    vector<Input> relayCells() const;
    vector<Input> interNeurons() const;
    vector<Input> corticalNeurons() const;

protected:
    Stimuli *m_stim;

    mat m_response, m_responseComplex;
    mat m_impulseResponse, m_impulseResponseComplex;

    vec m_mesh;
    vec3 m_domain;


    vector<Input> m_ganglionCells;
    vector<Input> m_relayCells;
    vector<Input> m_interNeurons;
    vector<Input> m_corticalNeurons;




};

#endif // NEURON_H
