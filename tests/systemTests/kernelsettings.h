#ifndef KERNELSETTINGS_H
#define KERNELSETTINGS_H

#include "kernels/temporalKernels/temporalkernel.h"
#include "kernels/spatialKernels/spatialkernel.h"

using namespace lgnSimulator;

class KernelSettings
{
public:
    KernelSettings();
    ~KernelSettings();

    static vector<SpatialKernel*> spatialKernelVector();
    static vector<TemporalKernel*> temporalKernelVector();
};

#endif // KERNELSETTINGS_H
