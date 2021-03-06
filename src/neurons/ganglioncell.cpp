#include "ganglioncell.h"

/*!
  \class lgnSimulator::GanglionCell
  \inmodule lgnSimulator
  \ingroup lgnSimulator-neurons
  \brief The GanglionCell class.
 */

using namespace lgnSimulator;


GanglionCell::GanglionCell(Integrator* const integrator,
                           const Kernel &kernel,
                           const double backgroundResponse,
                           StaticNonlinearity *staticNonlinearity)
    : Neuron(integrator, backgroundResponse, "ganglion", staticNonlinearity)
    , m_kernel(kernel)
{
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
                {m_spatialFreqs[i], m_spatialFreqs[j]}, m_temporalFreqs[k]);
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
    return m_kernel.spatiotemporal(rVec,t);
}

complex<double> GanglionCell::impulseResponseFourierTransformAtFrequency(vec2 kVec,
                                                                         double w)
{
    return m_kernel.fourierTransform(kVec,w);
}

