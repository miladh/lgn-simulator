#ifndef NEURON_H
#define NEURON_H

#include <armadillo>
#include <libconfig.h++>
#include <fftw3.h>


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


    // Compute functions
    void computeResponse(double t);
    void computeResponseComplex(double w);
    void computeImpulseResponse(double t);
    void computeImpulseResponseComplex(double w);

    // Virtual functions
    virtual double impulseResponseComplex(vec2 kVec, double w) = 0;

    // Add cell functions
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
     cx_mat complexResponse ;

    vec m_spatialMesh;
    vec m_freqMesh;
    vec3 m_domain;


    vector<Input> m_ganglionCells;
    vector<Input> m_relayCells;
    vector<Input> m_interNeurons;
    vector<Input> m_corticalNeurons;

    const std::complex<double> m_i = {0,1};

};

#endif // NEURON_H
