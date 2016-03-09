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



complex<double> RelayCell::impulseResponseFourierTransformAtFrequency(int idx,
                                                                      int jdx,
                                                                      int kdx)
{
    vec2 kVec= {m_spatialFreqs[idx], m_spatialFreqs[jdx]};
    double w = m_temporalFreqs[kdx];


    complex<double> ganglionFF = 0;
    complex<double> interneuronFF = 0;
    complex<double> interneuronFB = 0;
    complex<double> corticalFB = 0;

    //Feedforward ganglion input
    for (const Input g : m_ganglionCells){
        Neuron* const ganglionCell = g.neuron;
        ganglionFF += g.kernel.fourierTransform(kVec,w)
                * ganglionCell->impulseResponseFourierTransform()(idx,jdx,kdx);
    }


    //Interneuron input
    for (const Input i : m_interNeurons){
        Interneuron* const interneuron  = dynamic_cast<Interneuron* const>(i.neuron);
        complex<double> Kri = i.kernel.fourierTransform(kVec,w);


        //Feedforward term
        for (const Input g : interneuron->ganglionCells()){
            Neuron* const ganglionCell = g.neuron;
            interneuronFF += g.kernel.fourierTransform(kVec,w)
                    * ganglionCell->impulseResponseFourierTransform()(idx,jdx,kdx);
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


    //Feedback cortical input
    complex<double> Kcr = 0.0;
    complex<double> Krc = 0.0;
    for (const Input c : m_corticalNeurons){
        CorticalCell* const corticalCell  =
                dynamic_cast<CorticalCell* const>(c.neuron);

        Krc = c.kernel.fourierTransform(kVec,w);
        const Kernel* kernel =  corticalCell->relayInputKernel();
        Kcr = kernel->fourierTransform(kVec,w);
        corticalFB += Krc * Kcr;
    }


    if((1 - interneuronFB.real()  - corticalFB.real() )== 0){
        cerr << "interneuronFB: " << interneuronFB << " corticalFB: " << corticalFB << endl;
        throw overflow_error("Divide by zero exception in relay feedback contribution");
    }


    //    cout << idx << "  " << jdx << "   " << kdx << endl;
    //    cout << ganglionFF << "  " << interneuronFF
    //         << "   " << interneuronFB << "   " << corticalFB << endl;
    //    cout << endl;

    complex<double> Wr = (ganglionFF + interneuronFF)
            /(complex<double>(1.0, 0.0) - interneuronFB - corticalFB);

    return Wr;
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

void RelayCell::addCorticalNeuron(Neuron* const neuron, const Kernel &kernel)
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





































