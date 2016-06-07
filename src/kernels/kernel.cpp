#include "kernel.h"

using namespace lgnSimulator;
Kernel::Kernel(double weight):
    m_weight(weight)

{

}

double Kernel::weight() const
{
    return m_weight;
}
