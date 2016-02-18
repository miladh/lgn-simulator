#ifndef GANGLIONCELL_H
#define GANGLIONCELL_H

#include "neuron.h"

#include <armadillo>

using namespace arma;
namespace lgnSimulator {
class GanglionCell : public Neuron
{
public:
    GanglionCell(Integrator *integrator,
                 Kernel *kernel,
                 double backgroundResponse= 0,
                 StaticNonlinearity *staticNonlinearity = nullptr);
    ~GanglionCell();

    virtual void computeImpulseResponse();
    virtual void computeImpulseResponseFourierTransform();

private:
    Kernel *m_kernel;

    double impulseResponseValueAtPoint(vec2 rVec, double t);
    complex<double> impulseResponseFourierTransformAtFrequency(vec2 kVec, double w);

};
}
#endif // GANGLIONCELL_H
