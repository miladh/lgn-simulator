#include "ganglioncell.h"

using namespace lgnSimulator;


GanglionCell::GanglionCell(Integrator *integrator,
                           SpatialKernel *spatialKernel,
                           TemporalKernel *temporalKernel,
                           double backgroundResponse,
                           StaticNonlinearity *staticNonlinearity)
    : Neuron(integrator, backgroundResponse, staticNonlinearity)
    , m_spatialKernel(spatialKernel)
    , m_temporalKernel(temporalKernel)
{
    m_cellType = "ganglion";
}

GanglionCell::~GanglionCell()
{

}

void GanglionCell::computeImpulseResponseFourierTransform()
{
    impulseResponseFourierTransformComputed = true;

    for(int k = 0; k < int(m_impulseResponseFT.n_slices); k++){
        for(int i = 0; i < int(m_impulseResponseFT.n_rows); i++){
            for(int j = 0; j < int(m_impulseResponseFT.n_cols); j++){
                m_impulseResponseFT(i,j,k) =
                        impulseResponseFourierTransformAtFrequency(
                {m_spatialFreqs[i], m_spatialFreqs[j]}, -m_temporalFreqs[k]);
            }
        }
    }
}


void GanglionCell::computeImpulseResponse()
{

    m_impulseResponse.set_size(m_impulseResponseFT.n_rows,
                               m_impulseResponseFT.n_rows,
                               m_impulseResponseFT.n_slices);

    for(int k = 0; k < int(m_impulseResponse.n_slices); k++){
        for(int i = 0; i < int(m_impulseResponse.n_rows); i++){
            for(int j = 0; j < int(m_impulseResponse.n_cols); j++){
                m_impulseResponse(i, j, k) =
                        impulseResponseValueAtPoint(
                {m_spatialVec[i], m_spatialVec[j]}, m_timeVec[k]);
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

