#ifndef GANGLIONCELL_H
#define GANGLIONCELL_H

#include "neuron.h"
#include "spatialKernels/spatialkernel.h"
#include "temporalKernels/temporalkernel.h"

class GanglionCell : public Neuron
{
public:
    GanglionCell(const Config *cfg, Stimuli *stim,
                 SpatialKernel *spatialKernel,
                 TemporalKernel *temporalKernel);
    ~GanglionCell();

private:
    SpatialKernel *m_spatialKernel;
    TemporalKernel *m_temporalKernel;


    double impulseResponse(vec2 rVec, double t);

    // Neuron interface
    double impulseResponseComplex(vec2 kVec, double w);
};

#endif // GANGLIONCELL_H
