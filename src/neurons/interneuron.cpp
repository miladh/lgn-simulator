#include "interneuron.h"


using namespace lgnSimulator;


Interneuron::Interneuron(const Integrator &integrator, double backgroundResponse)
    : Neuron(integrator, backgroundResponse, "interneuron")
{
}

Interneuron::~Interneuron()
{

}

void Interneuron::computeImpulseResponseFourierTransform()
{
    computeNeededcubes();
    impulseResponseFourierTransformComputed = true;

    for(int k = 0; k < int(m_impulseResponseFT.n_slices); k++){
        double w = m_temporalFreqs[k];

        for(int i = 0; i < int(m_impulseResponseFT.n_rows); i++){
            for(int j = 0; j < int(m_impulseResponseFT.n_cols); j++){
                vec2 kVec= {m_spatialFreqs[i], m_spatialFreqs[j]};

                for (const Input g : m_ganglionCells){
                    Neuron* const ganglionCell = g.neuron;
                    m_impulseResponseFT(i,j,k)
                            += g.kernel.fourierTransform(kVec,w)
                            * ganglionCell->
                            impulseResponseFourierTransform()(i,j,k);
                }

                for (const Input c : m_corticalNeurons){
                    Neuron* const corticalCell = c.neuron;
                    m_impulseResponseFT(i,j,k)
                            += c.kernel.fourierTransform(kVec,w)
                            * corticalCell->
                            impulseResponseFourierTransform()(i,j,k);
                }
            }
        }
    }
}

void Interneuron::computeNeededcubes() const
{
    for (const Input g : m_ganglionCells){
        Neuron* const ganglionCell = g.neuron;
        if(!ganglionCell->isImpulseResponseFourierTransformComputed()){
            ganglionCell->computeImpulseResponseFourierTransform();
        }
    }

    for (const Input c : m_corticalNeurons){
        Neuron* const corticalCell = c.neuron;
        if(!corticalCell->isImpulseResponseFourierTransformComputed()){
            corticalCell->computeImpulseResponseFourierTransform();
        }
    }
}



void Interneuron::addGanglionCell(Neuron* const neuron, const Kernel &kernel)
{
    if (neuron->type() == "ganglion") {
        m_ganglionCells.emplace_back(Input{neuron, kernel});
    }else{
        throw overflow_error("wrong cell type in addGanglionCell(): "
                             + neuron->type());
    }

}


void Interneuron::addCorticalCell(Neuron* const neuron, const Kernel &kernel)
{
    if (neuron->type() == "cortical") {
        m_corticalNeurons.emplace_back(Input{neuron, kernel});
    }else{throw overflow_error("wrong cell type in addCorticalNeuron(): "
                             + neuron->type());
    }
}


vector<Neuron::Input> Interneuron::ganglionCells() const
{
    return m_ganglionCells;
}

vector<Neuron::Input> Interneuron::corticalNeurons() const
{
    return m_corticalNeurons;
}








