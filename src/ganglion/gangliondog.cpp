#include "gangliondog.h"

GanglionDOG::GanglionDOG(DOG *dog)
    : Ganglion()
    , m_dog(dog)
{

}

GanglionDOG::~GanglionDOG()
{

}

void GanglionDOG::setSpatialImpulseResponseReal(vec2 r)
{
    m_spatialImpulseResponseReal = m_dog->real(r);
}

void GanglionDOG::setSpatialImpulseResponseComplex(vec2 k)
{

    m_spatialImpulseResponseComplex = m_dog->complex(k);
}

void GanglionDOG::setTemporalImpulseResponseReal(double t)
{
    m_temporalImpulseResponseReal = 1.;
}

void GanglionDOG::setTemporalImpulseResponseComplex(double w)
{
    m_temporalImpulseResponseComplex = 1.0;
}

