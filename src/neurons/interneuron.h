#ifndef INTERNEURONS_H
#define INTERNEURONS_H

#include "neuron.h"

namespace lgnSimulator {
class Interneuron : public Neuron
{

public:
    Interneuron(const Integrator &integrator, double backgroundResponse= 0);
    ~Interneuron();

    void addGanglionCell(Neuron* const neuron, const Kernel &kernel);
    void addCorticalNeuron(Neuron* const neuron, const Kernel &kernel);

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
