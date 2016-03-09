#include "kernelsettings.h"

#include "kernels/separablekernel.h"
#include "kernels/spatialKernels/dog.h"
#include "kernels/spatialKernels/ellipticgaussian.h"

#include "kernels/temporalKernels/decayingexponential.h"
#include "kernels/temporalKernels/biphasic.h"
#include "kernels/temporalKernels/temporallyconstant.h"

KernelSettings::KernelSettings()
{

}

KernelSettings::~KernelSettings()
{
}

vector<SpatialKernel *> KernelSettings::spatialKernelVector()
{

    //DOG
    double A = 1.0;
    double a = 2.1;
    double B = 0.5;
    double b = 0.8;

    //Elliptic gauss
    double weightElliptic = 2.0;
    double angle = 45.0;
    double widthLong = 0.1;
    double widthNarrow = 1.0;

    vector<SpatialKernel*> spatialKernels = {new DOG(A,a,B,b),
                                             new EllipticGaussian
                                             (weightElliptic, angle, widthLong,
                                             widthNarrow)};

    return spatialKernels;

}

vector<TemporalKernel *> KernelSettings::temporalKernelVector()
{
    //Decaying exponential:
    double tau = 1.0;
    double delay = 0.0;


    //Damped oscillator
    double phaseDuration = 1.0;
    double weightDampedOsc = 2.0;

    vector<TemporalKernel*> temporalKernels = {new DecayingExponential(tau, delay),
                                              new Biphasic
                                              (phaseDuration, weightDampedOsc, delay)};

    return temporalKernels;
}

