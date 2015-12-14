#include "heavisidenonlinearity.h"

HeavisideNonlinearity::HeavisideNonlinearity()
{

}

HeavisideNonlinearity::~HeavisideNonlinearity()
{

}


double HeavisideNonlinearity::advance(const double u)
{
    return Functions::heaviside(u);
}
