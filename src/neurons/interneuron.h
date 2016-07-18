#ifndef INTERNEURONS_H
#define INTERNEURONS_H

#include "neuron.h"

namespace lgnSimulator {
class Interneuron : public Neuron
{

public:
    Interneuron(Integrator* const integrator, double backgroundResponse= 0);
    ~Interneuron();

    void addGanglionCell(Neuron* const neuron, const Kernel &kernel);
    void addCorticalCell(Neuron* const neuron, const Kernel &kernel);

    vector<Input> ganglionCells() const;
    vector<Input> corticalNeurons() const;

    // Neuron interface
    virtual void computeImpulseResponseFourierTransform();

private:
    vector<Input> m_ganglionCells;
    vector<Input> m_corticalNeurons;

    void computeNeededcubes() const;
};
}
#endif // INTERNEURONS_H
