#include "relaycell.h"

using namespace lgnSimulator;

/*!
 * \class RelayCell
 * \inmodule lgnSimulator
 * \ingroup lgnSimulator-neurons
 * \brief Relay cell class.
 *
 */

RelayCell::RelayCell(const Integrator &integrator, double backgroundResponse)
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

    for (const Input i : m_interNeurons){
        Interneuron* const interneuron  = dynamic_cast<Interneuron* const>(i.neuron);
        for (const Input g : interneuron->ganglionCells()){
            Neuron *ganglionCell = g.neuron;
            if(!ganglionCell->isImpulseResponseFourierTransformComputed()){
                ganglionCell->computeImpulseResponseFourierTransform();
            }
        }
    }
}



complex<double> RelayCell::impulseResponseFourierTransformAtFrequency(int kxi,
                                                                      int kyi,
                                                                      int wi) const
{
    cx_vec interneuron =  interneuronInput(kxi, kyi, wi);
    complex<double> ganglionFF = ganglionInput(kxi, kyi, wi);
    complex<double> interneuronFF = interneuron(0);
    complex<double> interneuronFB = interneuron(1);
    complex<double> corticalFB = corticalInput(kxi, kyi, wi);

    complex<double> fb = 1. - interneuronFB  - corticalFB;
    if(fabs(fb.real()) < 1.e-10 && fabs(fb.imag()) < 1.e-10){
        cerr << "interneuronFB: " << interneuronFB
             << " corticalFB: "   << corticalFB   << endl;
        throw overflow_error("Divide by zero exception in relay feedback contribution");
    }


    complex<double> Wr = (ganglionFF + interneuronFF)/fb;

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



cx_vec RelayCell::interneuronInput(int kxi, int kyi, int wi) const
{
    vec2 kVec= {m_spatialFreqs[kxi], m_spatialFreqs[kyi]};
    double w = m_temporalFreqs[wi];

    complex<double> interneuronFF = complex<double>(0, 0);;
    complex<double> interneuronFB = complex<double>(0, 0);;


    //Interneuron input
    for (const Input i : m_interNeurons){
        Interneuron* const interneuron  = dynamic_cast<Interneuron* const>(i.neuron);
        complex<double> Kri = i.kernel.fourierTransform(kVec,w);


        //Feedforward term
        for (const Input g : interneuron->ganglionCells()){
            Neuron* const ganglionCell = g.neuron;
            interneuronFF += g.kernel.fourierTransform(kVec,w)
                    * ganglionCell->impulseResponseFourierTransform()(kxi, kyi, wi);
        }
        interneuronFF*= Kri;

        //Feedback term
        for (const Input c : interneuron->corticalNeurons()){
            CorticalCell* const corticalCell  =
                    dynamic_cast<CorticalCell* const>(c.neuron);
            complex<double> Kic = c.kernel.fourierTransform(kVec,w);

            const Kernel* kernel =  corticalCell->relayInputKernel();
            complex<double> Kcr = kernel->fourierTransform(kVec,w);
            interneuronFB += Kri*Kic*Kcr;

        }

    }

    return cx_vec{interneuronFF, interneuronFB};
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


void RelayCell::addInterNeuron(Neuron* const neuron, const Kernel &kernel)
{
    if (neuron->type() == "interneuron") {
        m_interNeurons.emplace_back(Input{neuron, kernel});
    }else{
        throw overflow_error("wrong cell type in addInterNeuron(): " + neuron->type());
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


vector<Neuron::Input> RelayCell::interNeurons() const
{
    return m_interNeurons;
}

vector<Neuron::Input> RelayCell::corticalNeurons() const
{
    return m_corticalNeurons;
}












