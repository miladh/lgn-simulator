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

double GanglionCell::impulseResponse(vec2 rVec, double t)
{
    return m_spatialKernel->real(rVec) * m_temporalKernel->real(t);
}

double GanglionCell::impulseResponseComplex(vec2 kVec, double w)
{
    return m_spatialKernel->complex(kVec) * m_temporalKernel->complex(w);
}

