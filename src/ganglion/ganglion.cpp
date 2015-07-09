#include "ganglion.h"

Ganglion::Ganglion()
{

}

Ganglion::~Ganglion()
{

}

double Ganglion::impulseResponseReal()
{
    return m_spatialImpulseResponseReal + m_temporalImpulseResponseReal;
}

double Ganglion::impulseResponseComplex()
{
    return m_spatialImpulseResponseComplex + m_temporalImpulseResponseComplex;
}


double Ganglion::spatialImpulseResponseReal() const
{
    return m_spatialImpulseResponseReal;
}

double Ganglion::spatialImpulseResponseComplex() const
{
    return m_spatialImpulseResponseComplex;
}

double Ganglion::temporalImpulseResponseReal() const
{
    return m_temporalImpulseResponseReal;
}

double Ganglion::temporalImpulseResponseComplex() const
{
    return m_temporalImpulseResponseComplex;
}







