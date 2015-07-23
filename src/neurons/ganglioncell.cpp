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

void GanglionCell::computeImpulseResponse(double t)
{
    m_impulseResponse  = 0* m_impulseResponse;

    for(int i = 0; i < int(m_mesh.n_elem); i++){
        for(int j = 0; j < int(m_mesh.n_elem); j++){
            m_impulseResponse(i, j) = impulseResponse({m_mesh[i], m_mesh[j]}, t);
        }
    }

}


double GanglionCell::impulseResponse(vec2 rVec, double t)
{
    return m_spatialKernel->real(rVec) * m_temporalKernel->real(t);
}

double GanglionCell::impulseResponseComplex(vec2 kVec, double w)
{
    return m_spatialKernel->complex(kVec) * m_temporalKernel->complex(w);
}

