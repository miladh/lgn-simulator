#include "ganglioncell.h"

GanglionCell::GanglionCell(const Config *cfg, Stimuli *stim,
                           Integrator integrator,
                           SpatialKernel *spatialKernel,
                           TemporalKernel *temporalKernel)
    : Neuron(cfg, stim, integrator)
    , m_spatialKernel(spatialKernel)
    , m_temporalKernel(temporalKernel)
{
    m_cellType = "ganglion";

}

GanglionCell::~GanglionCell()
{

}

void GanglionCell::computeImpulseResponse()
{
    for(int k = 0; k < m_impulseResponse.n_slices; k++){
        for(int i = 0; i < m_impulseResponse.n_rows; i++){
            for(int j = 0; j < m_impulseResponse.n_cols; j++){
                m_impulseResponse(i, j, k) =
                        impulseResponse({m_coordinateVec[i]-0.5, m_coordinateVec[j]-0.5},
                                        timeVec[k]);
            }
        }
    }
}


double GanglionCell::impulseResponse(vec2 rVec, double t)
{
    return m_spatialKernel->real(rVec) * m_temporalKernel->real(t);
}

double GanglionCell::impulseResponseFT(vec2 kVec, double w)
{
    return m_spatialKernel->complex(kVec) * m_temporalKernel->complex(w);
}

