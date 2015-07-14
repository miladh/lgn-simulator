#include "relaycell.h"

RelayCell::RelayCell(const Config *cfg)
    : Neuron(cfg)
{

}

RelayCell::~RelayCell()
{

}


void RelayCell::computeResponse(double t)
{
    mat stim = 0*m_response;
    m_response = 0*m_response;
    m_impulseResponse  = 0* m_impulseResponse;

    double *w = new double [int(m_domain[2])];
    double *x = new double [int(m_domain[2])];
    gauleg(m_domain[0], m_domain[1], x, w, m_domain[2]);

    for(int i = 0; i < int(m_mesh.n_elem); i++){
        for(int j = 0; j < int(m_mesh.n_elem); j++){

            stim(i,j) = m_stim->real({m_mesh[i], m_mesh[j]}, t);
            for(int m = 0; m < int(m_domain[2]); m++){
                for(int n = 0; n < int(m_domain[2]); n++){

                    double Scomplex = m_stim->complex({x[m], x[n]}, m_stim->w());
                    double Gcomplex = impulseResponseComplex({x[m], x[n]}, m_stim->w());

                    double dGr = Gcomplex * w[m] * w[n] *
                            cos(m_mesh[i]*x[m]+ m_mesh[j]*x[n] - m_stim->w() * t);

                    m_impulseResponse(i, j) +=  dGr;
                    m_response(i,j) += dGr * Scomplex;
                }
            }
        }
    }
    m_stim->setReal(stim);
}

double RelayCell::impulseResponseComplex(vec2 kVec, double w)
{

    double G, I, C = 0;

    for (const Input g : m_ganglionCells){
        Neuron *ganglionCell = g.neuron;
        G += g.spatialKernel->coupling(kVec,w)
           * g.temporalKernel->coupling(kVec,w)
           * ganglionCell->impulseResponseComplex(kVec,w);
    }

    for (const Input i : m_interNeurons){
        Neuron *interneuron = i.neuron;
        double Kri = i.spatialKernel->coupling(kVec,w)
                * i.temporalKernel->coupling(kVec,w);

        for (const Input g : interneuron->ganglionCells()){
            Neuron *ganglionCell = g.neuron;
            I += g.spatialKernel->coupling(kVec,w)
               * g.temporalKernel->coupling(kVec,w)
               * ganglionCell->impulseResponseComplex(kVec,w);
        }
        I*= Kri;
    }


    for (const Input c : m_corticalNeurons){
        Neuron *corticalCell = c.neuron;
        double Krc = c.spatialKernel->coupling(kVec,w)
                * c.temporalKernel->coupling(kVec,w);

        for (const Input r : corticalCell->relayCells()){
            C += r.spatialKernel->coupling(kVec,w)
               * r.temporalKernel->coupling(kVec,w);
        }
        C*= Krc;
    }


    double Gr = (G + I)/(1 - C);
    return Gr;
}

