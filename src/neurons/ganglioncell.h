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

    // Neuron interface
    void computeResponse(double t);
    void computeResponseComplex(double w);
    double impulseResponse(vec2 rVec, double t);
    double impulseResponseComplex(vec2 kVec, double w);

private:
    SpatialKernel *m_spatialKernel;
    TemporalKernel *m_temporalKernel;

};

#endif // GANGLIONCELL_H
