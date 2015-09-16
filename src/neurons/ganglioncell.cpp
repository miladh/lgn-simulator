#include "ganglioncell.h"

GanglionCell::GanglionCell(const Config *cfg, Stimuli *stim,
                           SpatialKernel *spatialKernel,
                           TemporalKernel *temporalKernel)
    : Neuron(cfg, stim)
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
    for(int k = 0; k < m_nSteps; k++){
        for(int i = 0; i < m_nPoints; i++){
            for(int j = 0; j < m_nPoints; j++){
                m_impulseResponse(i, j, k) =
                        impulseResponse({m_spatialMesh[i], m_spatialMesh[j]},
                                        m_temporalMesh[k]);
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

