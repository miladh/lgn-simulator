#ifndef TEMPORALKERNELS_H
#define TEMPORALKERNELS_H

#include <armadillo>

#include "integrator.h"
#include "helper/specialfunctions.h"

using namespace std;
using namespace arma;
namespace lgnSimulator {
class TemporalKernel
{
public:
    TemporalKernel();
    ~TemporalKernel();

    virtual double temporal(double t) const = 0;
    virtual complex<double> fourierTransform(double w) const = 0;

};
}
#endif // TEMPORALKERNELS_H
