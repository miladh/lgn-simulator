#include "ganglioncell.h"

GanglionCell::GanglionCell(const Config *cfg, Stimuli *stim, SpatialKernel *spatialKernel, TemporalKernel *temporalKernel)
    : Neuron(cfg, stim)
    , m_spatialKernel(spatialKernel)
    , m_temporalKernel(temporalKernel)
{
    m_cellType = "ganglion";
}

GanglionCell::~GanglionCell()
{

}

void GanglionCell::computeResponse(double t)
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

void GanglionCell::computeResponseComplex(double w)
{
    for(int i = 0; i < int(m_mesh.n_elem); i++){
        for(int j = 0; j < int(m_mesh.n_elem); j++){

            double Gcomplex = impulseResponseComplex({m_mesh[i], m_mesh[j]}, w);
            double Scomplex = m_stim->complex({m_mesh[i], m_mesh[j]}, w);

            m_impulseResponseComplex(i,j) = Gcomplex;
            m_responseComplex(i,j) = Gcomplex*Scomplex;
        }
    }
}

double GanglionCell::impulseResponseComplex(vec2 kVec, double w)
{
    return m_spatialKernel->complex(kVec) * m_temporalKernel->complex(w);
}

