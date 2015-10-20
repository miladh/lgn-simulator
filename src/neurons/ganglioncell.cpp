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
    for(int k = 0; k < int(m_impulseResponse.n_slices); k++){
        for(int i = 0; i < int(m_impulseResponse.n_rows); i++){
            for(int j = 0; j < int(m_impulseResponse.n_cols); j++){
                m_impulseResponse(i, j, k) =
                        impulseResponse({m_coordinateVec[i], m_coordinateVec[j]},
                                        timeVec[k]);
            }
        }
    }
}


double GanglionCell::impulseResponse(vec2 rVec, double t)
{
    return m_spatialKernel->spatial(rVec) * m_temporalKernel->temporal(t);
}

double GanglionCell::impulseResponseFT(vec2 kVec, double w)
{
    return m_spatialKernel->fourierTransform(kVec) * m_temporalKernel->fourierTransform(w);
}

