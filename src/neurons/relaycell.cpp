#include "relaycell.h"

RelayCell::RelayCell(Integrator *integrator)
    : Neuron(integrator)

{
    m_cellType = "relay";
}

RelayCell::~RelayCell()
{

}


complex<double> RelayCell::impulseResponseFourierTransformAtFrequency(vec2 kVec, double w)
{

    complex<double> G = 0;
    complex<double> Iff = 0;
    complex<double> Ifb = 0;
    complex<double> C = 0;

    for (const Input g : m_ganglionCells){
        Neuron *ganglionCell = g.neuron;
        G += g.spatialKernel->fourierTransform(kVec)
                * g.temporalKernel->fourierTransform(w)
                * ganglionCell->impulseResponseFourierTransformAtFrequency(kVec,w);
    }



    for (const Input i : m_interNeurons){
        Neuron *interneuron = i.neuron;
        complex<double> Kri = i.spatialKernel->fourierTransform(kVec)
                * i.temporalKernel->fourierTransform(w);


        //Feedforward term
        for (const Input g : interneuron->ganglionCells()){
            Neuron *ganglionCell = g.neuron;
            Iff += g.spatialKernel->fourierTransform(kVec)
                    * g.temporalKernel->fourierTransform(w)
                    * ganglionCell->impulseResponseFourierTransformAtFrequency(kVec,w);
        }
        Iff*= Kri;

        //Feedback term
        for (const Input c : interneuron->corticalNeurons()){
            Neuron *corticalCell = c.neuron;
            complex<double> Kic = c.spatialKernel->fourierTransform(kVec)
                    * c.temporalKernel->fourierTransform(w);

            // NOTE: ONLY ONE RELAY CELL!!!
            for (const Input r : corticalCell->relayCells()){
                complex<double> Kcr = r.spatialKernel->fourierTransform(kVec)
                        * r.temporalKernel->fourierTransform(w);
                Ifb += Kri*Kic*Kcr;
            }

        }

    }



    for (const Input c : m_corticalNeurons){
        Neuron *corticalCell = c.neuron;
        complex<double> Krc = c.spatialKernel->fourierTransform(kVec)
                * c.temporalKernel->fourierTransform(w);

        for (const Input r : corticalCell->relayCells()){
            C += r.spatialKernel->fourierTransform(kVec)
                    * r.temporalKernel->fourierTransform(w);
        }
        C*= Krc;
    }

    if((1 - Ifb.real()  - C.real() )== 0){
        cout << "Ifb: " << Ifb << " C: " << C << endl;
        throw overflow_error("Divide by zero exception in relay feedback contribution");
    }

    complex<double> Gr = (G + Iff)/(complex<double>(1.0, 0.0) - Ifb - C);

    return Gr;
}

