#include "kernel.h"


/*!
 \class lgnSimulator::Kernel
 \inmodule lgnSimulator
 \brief The Kernel class.
 */


using namespace lgnSimulator;
Kernel::Kernel(double weight):
    m_weight(weight)

{

}

double Kernel::weight() const
{
    return m_weight;
}
