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

    Neuron * ganglion = m_feedforwardInputs[0].neuron;
    for(int i = 0; i < int(m_mesh.n_elem); i++){
        for(int j = 0; j < int(m_mesh.n_elem); j++){

            stim(i,j) = m_stim->real({m_mesh[i], m_mesh[j]}, t);
            for(int m = 0; m < int(m_domain[2]); m++){
                for(int n = 0; n < int(m_domain[2]); n++){

                    double T = transferFunctionComplex({x[m], x[n]}, m_stim->w());
//                    double Gg = ganglion->impulseResponseComplex({x[m], x[n]}, m_stim->w());
                    double Gg = 1;
                    double s = m_stim->complex({x[m], x[n]}, m_stim->w());

                    double dGr = T * Gg * w[m] * w[n] *
                            cos(m_mesh[i]*x[m]+ m_mesh[j]*x[n] - m_stim->w() * t);
                    m_impulseResponse(i, j) +=  dGr;
                    m_response(i,j) += dGr * s;
                }
            }
        }
    }
    m_stim->setReal(stim);
}

double RelayCell::transferFunctionComplex(vec2 kVec, double w)
{

    double ff = 0;
    for (const Input inp : m_feedforwardInputs){
        ff += inp.spatialKernel->coupling(kVec,w) * inp.temporalKernel->coupling(kVec,w);
    }

    double fb = 0;
    for (const Input inp : m_feedbackInputs){
        fb += inp.spatialKernel->coupling(kVec,w) * inp.temporalKernel->coupling(kVec,w);
    }
    return ff / (1 - fb);
}

