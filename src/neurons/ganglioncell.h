#ifndef GANGLIONCELL_H
#define GANGLIONCELL_H

#include "neuron.h"
#include "spatialKernels/spatialkernel.h"
#include "temporalKernels/temporalkernel.h"

#include <armadillo>

using namespace arma;

class GanglionCell : public Neuron
{
public:
    GanglionCell(Stimuli *stim, Integrator integrator,
                 SpatialKernel *spatialKernel,
                 TemporalKernel *temporalKernel);
    ~GanglionCell();

    void computeImpulseResponse();

private:
    SpatialKernel *m_spatialKernel;
    TemporalKernel *m_temporalKernel;


    double impulseResponse(vec2 rVec, double t);

    // Neuron interface
    double impulseResponseFT(vec2 kVec, double w);
};

#endif // GANGLIONCELL_H
