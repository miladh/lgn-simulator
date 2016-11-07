#include "relaycell.h"

using namespace lgnSimulator;

/*!
  \class lgnSimulator::RelayCell
  \inmodule lgnSimulator
  \ingroup lgnSimulator-neurons
  \brief The RelayCell class.
 */

RelayCell::RelayCell(Integrator* const integrator, double backgroundResponse)
    : Neuron(integrator, backgroundResponse,  "relay")

{
}

RelayCell::~RelayCell()
{
}

void RelayCell::computeImpulseResponseFourierTransform()
{
    computeNeededcubes();
    impulseResponseFourierTransformComputed = true;

    for(int k = 0; k < int(m_impulseResponseFT.n_slices); k++){
        for(int i = 0; i < int(m_impulseResponseFT.n_rows); i++){
            for(int j = 0; j < int(m_impulseResponseFT.n_cols); j++){
                m_impulseResponseFT(i,j,k) =
                        impulseResponseFourierTransformAtFrequency(i, j, k);
            }
        }
    }
}

void RelayCell::computeNeededcubes() const
{
    for (const Input g : m_ganglionCells){
        Neuron *ganglionCell = g.neuron;
        if(!ganglionCell->isImpulseResponseFourierTransformComputed()){
            ganglionCell->computeImpulseResponseFourierTransform();
        }
    }
}


complex<double> RelayCell::impulseResponseFourierTransformAtFrequency(int kxi,
                                                                      int kyi,
                                                                      int wi) const
{
    complex<double> ganglionFF = ganglionInput(kxi, kyi, wi);
    complex<double> corticalFB = corticalInput(kxi, kyi, wi);

    complex<double> fb = 1. -  corticalFB;
    if(fabs(fb.real()) < 1.e-10 && fabs(fb.imag()) < 1.e-10){
        cerr << " corticalFB: "   << corticalFB   << endl;
        throw overflow_error("Divide by zero exception in relay feedback contribution");
    }


    complex<double> Wr = (ganglionFF)/fb;

    return Wr;
}



complex<double> RelayCell::ganglionInput(int kxi, int kyi, int wi) const
{
    vec2 kVec= {m_spatialFreqs[kxi], m_spatialFreqs[kyi]};
    double w = m_temporalFreqs[wi];

    complex<double> ganglionFF = complex<double>(0, 0);

    //Feedforward ganglion input
    for (const Input g : m_ganglionCells){
        Neuron* const ganglionCell = g.neuron;
        ganglionFF += g.kernel.fourierTransform(kVec,w)
                * ganglionCell->impulseResponseFourierTransform()(kxi, kyi, wi);
    }

    return ganglionFF;
}

complex<double> RelayCell::corticalInput(int kxi, int kyi, int wi)const
{
    vec2 kVec= {m_spatialFreqs[kxi], m_spatialFreqs[kyi]};
    double w = m_temporalFreqs[wi];

    complex<double> corticalFB = complex<double>(0, 0);

    //Feedback cortical input
    for (const Input c : m_corticalNeurons){
        CorticalCell* const corticalCell  =
                dynamic_cast<CorticalCell* const>(c.neuron);

        complex<double> Krc = c.kernel.fourierTransform(kVec,w);
        const Kernel* kernel =  corticalCell->relayInputKernel();
        complex<double> Kcr  = kernel->fourierTransform(kVec,w);

        corticalFB += Krc * Kcr;

    }

    return corticalFB;

}



void RelayCell::addGanglionCell(Neuron* const neuron, const Kernel &kernel)
{
    if (neuron->type() == "ganglion") {
        m_ganglionCells.emplace_back(Input{neuron, kernel});
    }else{
        throw overflow_error("wrong cell type in addGanglionCell(): " + neuron->type());
    }

}


void RelayCell::addCorticalCell(Neuron* const neuron, const Kernel &kernel)
{
    if (neuron->type() == "cortical") {
        m_corticalNeurons.emplace_back(Input{neuron, kernel});
    }else{
        throw overflow_error("wrong cell type in addCorticalNeuron(): " + neuron->type());
    }
}


vector<Neuron::Input> RelayCell::ganglionCells() const
{
    return m_ganglionCells;
}

vector<Neuron::Input> RelayCell::corticalNeurons() const
{
    return m_corticalNeurons;
}












