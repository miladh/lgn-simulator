#ifndef NEURON_H
#define NEURON_H


#include "../stimuli/stimuli.h"
#include "integrator.h"
#include "../temporalKernels/temporalkernel.h"
#include "../spatialKernels/spatialkernel.h"


using namespace arma;
using namespace std;

class Neuron
{
public:
    Neuron(Integrator *integrator);
    ~Neuron();

    struct Input {
        Neuron *neuron;
        SpatialKernel *spatialKernel;
        TemporalKernel *temporalKernel;
    };


    // Compute functions
    void computeResponse(Stimulus *stimulus);
    void computeImpulseResponse();


    // Virtual functions
    virtual complex<double> impulseResponseFourierTransformAtFrequency(vec2 kVec, double w) = 0;


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
    cube response() const;
    cube impulseResponse() const;
    cx_cube responseFT() const;
    cx_cube impulseResponseFourierTransformAtFrequency() const;


    vector<Input> ganglionCells() const;
    vector<Input> relayCells() const;
    vector<Input> interNeurons() const;
    vector<Input> corticalNeurons() const;

    string cellType() const;

private:
    Integrator* m_integrator;

protected:
    string m_cellType;

    cube m_response;
    cube m_impulseResponse;
    cx_cube m_responseFT;
    cx_cube m_impulseResponseFT;

    vec m_coordinateVec;
    vec m_spatialFreqs;

    vec timeVec;
    vec m_temporalFreqs;


    vector<Input> m_ganglionCells;
    vector<Input> m_relayCells;
    vector<Input> m_interNeurons;
    vector<Input> m_corticalNeurons;

    void computeImpulseResponseFourierTransform();
};

#endif // NEURON_H
