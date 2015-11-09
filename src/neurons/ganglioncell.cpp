#include "ganglioncell.h"

GanglionCell::GanglionCell(Integrator *integrator,
                           SpatialKernel *spatialKernel,
                           TemporalKernel *temporalKernel)
    : Neuron(integrator)
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
    computeImpulseResponseFourierTransform();
    for(int k = 0; k < int(m_impulseResponse.n_slices); k++){
        for(int i = 0; i < int(m_impulseResponse.n_rows); i++){
            for(int j = 0; j < int(m_impulseResponse.n_cols); j++){
                m_impulseResponse(i, j, k) =
                        impulseResponseValueAtPoint(
                {m_coordinateVec[i], m_coordinateVec[j]}, timeVec[k]);
            }
        }
    }
}


double GanglionCell::impulseResponseValueAtPoint(vec2 rVec, double t)
{
    return m_spatialKernel->spatial(rVec) * m_temporalKernel->temporal(t);
}

complex<double> GanglionCell::impulseResponseFourierTransformAtFrequency(vec2 kVec, double w)
{
    return m_spatialKernel->fourierTransform(kVec)
            * m_temporalKernel->fourierTransform(w);
}

