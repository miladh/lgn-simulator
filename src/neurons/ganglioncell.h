#ifndef GANGLIONCELL_H
#define GANGLIONCELL_H

#include "neuron.h"

class GanglionCell : public Neuron
{
public:
    GanglionCell(const Config *cfg, Stimuli *stim);
    ~GanglionCell();

    // Neuron interface
    void computeResponse(double t);
    void computeResponseComplex(double w);
    double impulseResponseComplex(vec2 kVec, double w);
};

#endif // GANGLIONCELL_H
