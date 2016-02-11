#include "relaycell.h"

using namespace lgnSimulator;

/*!
 * \class RelayCell
 * \inmodule lgnSimulator
 * \ingroup lgnSimulator-neurons
 * \brief Relay cell class.
 *
 */

RelayCell::RelayCell(Integrator *integrator)
    : Neuron(integrator)

{
    m_cellType = "relay";
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


void RelayCell::computeNeededcubes()
{
    for (const Input g : m_ganglionCells){
        Neuron *ganglionCell = g.neuron;
        if(!ganglionCell->isImpulseResponseFourierTransformComputed()){
            ganglionCell->computeImpulseResponseFourierTransform();
        }
    }

    for (const Input i : m_interNeurons){
        Neuron *interneuron = i.neuron;
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
    double w = -m_temporalFreqs[kdx];


    complex<double> ganglionFF = 0;
    complex<double> interneuronFF = 0;
    complex<double> interneuronFB = 0;
    complex<double> corticalFB = 0;

    //Feedforward ganglion input
    for (const Input g : m_ganglionCells){
        Neuron *ganglionCell = g.neuron;
        ganglionFF += g.spatialKernel->fourierTransform(kVec)
                * g.temporalKernel->fourierTransform(w)
                * ganglionCell->impulseResponseFourierTransform()(idx,jdx,kdx);
    }



    //Interneuron input
    for (const Input i : m_interNeurons){
        Neuron *interneuron = i.neuron;
        complex<double> Kri = i.spatialKernel->fourierTransform(kVec)
                * i.temporalKernel->fourierTransform(w);


        //Feedforward term
        for (const Input g : interneuron->ganglionCells()){
            Neuron *ganglionCell = g.neuron;
            interneuronFF += g.spatialKernel->fourierTransform(kVec)
                    * g.temporalKernel->fourierTransform(w)
                    * ganglionCell->impulseResponseFourierTransform()(idx,jdx,kdx);
        }
        interneuronFF*= Kri;

        //Feedback term
        for (const Input c : interneuron->corticalNeurons()){
            Neuron *corticalCell = c.neuron;
            complex<double> Kic = c.spatialKernel->fourierTransform(kVec)
                    * c.temporalKernel->fourierTransform(w);

            // NOTE: ONLY ONE RELAY CELL!!!
            for (const Input r : corticalCell->relayCells()){
                complex<double> Kcr = r.spatialKernel->fourierTransform(kVec)
                        * r.temporalKernel->fourierTransform(w);
                interneuronFB += Kri*Kic*Kcr;
            }

        }

    }


    //Feedback cortical input
    complex<double> Kcr = 0.0;
    complex<double> Krc = 0.0;
    for (const Input c : m_corticalNeurons){
        Neuron *corticalCell = c.neuron;
        Krc = c.spatialKernel->fourierTransform(kVec)
                * c.temporalKernel->fourierTransform(w);

        // NOTE: ONLY ONE RELAY CELL!!!
        for (const Input r : corticalCell->relayCells()){
            Kcr = r.spatialKernel->fourierTransform(kVec)
                    * r.temporalKernel->fourierTransform(w);
        }
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

