#ifndef STATICNONLINEARITY_H
#define STATICNONLINEARITY_H

#include <armadillo>
#include <yaml-cpp/yaml.h>

using namespace arma;

namespace lgnSimulator {
class StaticNonlinearity
{
public:
    StaticNonlinearity();
    ~StaticNonlinearity();

   void applyStaticNonlinearity(cube* const linearFilter);

protected:
    virtual double advance(const double u) = 0;
};
}
#endif // STATICNONLINEARITY_H
