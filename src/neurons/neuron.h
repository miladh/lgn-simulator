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
    Neuron(const Config *cfg, Stimuli *stim);
    ~Neuron();

    struct Input {
        Neuron *neuron;
        SpatialKernel *spatialKernel;
        TemporalKernel *temporalKernel;
    };


    void addGanglionCell(Neuron *neuron,
                         SpatialKernel *sKernel,
                         TemporalKernel *tKernel);

    void addRelayCell(Neuron *neuron,
                      SpatialKernel *sKernel,
                      TemporalKernel *tKernel);

    void addInterNeuron(Neuron *neuron,
                        SpatialKernel *sKernel,
                        TemporalKernel *tKernel);
    void addCorticalNeuron(Neuron *neuron,
                           SpatialKernel *sKernel,
                           TemporalKernel *tKernel);



    // Virtual functions
    virtual void computeResponse(double t) = 0;
    virtual void computeResponseComplex(double w) = 0;
    virtual double impulseResponseComplex(vec2 kVec, double w) = 0;




    // Getter member functions
    mat response() const;
    mat impulseResponse() const;
    mat responseComplex() const;
    mat impulseResponseComplex() const;


    vector<Input> ganglionCells() const;
    vector<Input> relayCells() const;
    vector<Input> interNeurons() const;
    vector<Input> corticalNeurons() const;

    string cellType() const;

protected:
    string m_cellType;
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
