#ifndef NEURON_H
#define NEURON_H


#include <libconfig.h++>
#include <fftw3.h>

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
    void computeResponse();
    void computeImpulseResponse();

    // Virtual functions
    virtual double impulseResponseFT(vec2 kVec, double w) = 0;

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
    cx_cube impulseResponseFT() const;


    vector<Input> ganglionCells() const;
    vector<Input> relayCells() const;
    vector<Input> interNeurons() const;
    vector<Input> corticalNeurons() const;

    string cellType() const;

private:
        void computeImpulseResponseFT();

protected:
    int m_nPoints = 0;
    int m_nSteps = 0;

    string m_cellType;
    Stimuli *m_stim;

    cube m_response;
    cube m_impulseResponse;
    cx_cube m_responseFT;
    cx_cube m_impulseResponseFT;

    vec m_spatialMesh;
    vec m_spatialFreqs;

    vec m_temporalMesh;
    vec m_temporalFreqs;


    vector<Input> m_ganglionCells;
    vector<Input> m_relayCells;
    vector<Input> m_interNeurons;
    vector<Input> m_corticalNeurons;

    const std::complex<double> m_i = {0,1};


};

#endif // NEURON_H
