#ifndef RELAYCELL_H
#define RELAYCELL_H

#include "neuron.h"
#include "corticalcell.h"

namespace lgnSimulator {

class CorticalCell;

class RelayCell : public Neuron
{
public:
    RelayCell(Integrator* const integrator, double backgroundResponse= 0);
    ~RelayCell();

    void addGanglionCell(Neuron* const neuron, const Kernel &kernel);
    void addCorticalCell(Neuron* const neuron, const Kernel &kernel);

    vector<Input> ganglionCells() const;
    vector<Input> corticalNeurons() const;

    // Neuron interface
    virtual void computeImpulseResponseFourierTransform();

private:
    vector<Input> m_ganglionCells;
    vector<Input> m_corticalNeurons;

    complex<double> ganglionInput(int kxi,int kyi, int wi) const;
    complex<double> corticalInput(int kxi,int kyi, int wi) const;
    complex<double> impulseResponseFourierTransformAtFrequency(int kxi,
                                                               int kyi,
                                                               int wi) const;

    void computeNeededcubes() const;
};
}
#endif // RELAYCELL_H
