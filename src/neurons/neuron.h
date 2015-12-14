#ifndef NEURON_H
#define NEURON_H


#include "../stimuli/stimuli.h"
#include "integrator.h"
#include "../temporalKernels/temporalkernel.h"
#include "../spatialKernels/spatialkernel.h"
#include "../staticNonlinearity/staticnonlinearity.h"


using namespace arma;
using namespace std;

class Neuron
{
public:
    Neuron(Integrator *integrator,
           StaticNonlinearity *staticNonlinearity = nullptr);
    ~Neuron();

    struct Input {
        Neuron *neuron;
        SpatialKernel *spatialKernel;
        TemporalKernel *temporalKernel;
    };


    // Compute functions
    void computeResponse(Stimulus *stimulus);
    void computeImpulseResponse();

    //Virtual function:
    virtual void computeImpulseResponseFourierTransform() = 0;


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
    const cube &response() const;
    const cube &impulseResponse() const;
    const cx_cube &responseFT() const;
    const cx_cube &impulseResponseFourierTransform() const;

    vector<Input> ganglionCells() const;
    vector<Input> relayCells() const;
    vector<Input> interNeurons() const;
    vector<Input> corticalNeurons() const;

    string cellType() const;


    void clearResponse();
    void clearImpulseResponse();

    bool isImpulseResponseFourierTransformComputed() const;

private:
    Integrator* m_integrator;
    StaticNonlinearity *m_staticNonlinearity;



protected:
    bool impulseResponseFourierTransformComputed  = false;
    string m_cellType;

    cube m_response;
    cube m_impulseResponse;
    cx_cube m_responseFT;
    cx_cube m_impulseResponseFT;

    vec m_coordinateVec;
    vec m_spatialFreqs;

    vec m_timeVec;
    vec m_temporalFreqs;


    vector<Input> m_ganglionCells;
    vector<Input> m_relayCells;
    vector<Input> m_interNeurons;
    vector<Input> m_corticalNeurons;

};

#endif // NEURON_H
