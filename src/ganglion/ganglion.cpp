#include "ganglion.h"

Ganglion::Ganglion(const Config *cfg)
{

}

Ganglion::~Ganglion()
{

}

double Ganglion::impulseResponse(vec2 r, double t)
{
    setSpatialImpulseResponse(r);
    setTemporalImpulseResponse(t);
    return m_spatialImpulseResponse + m_temporalImpulseResponse;
}

double Ganglion::impulseResponseComplex(vec2 k, double w)
{
    setSpatialImpulseResponseComplex(k);
    setTemporalImpulseResponseComplex(w);
    return m_spatialImpulseResponseComplex + m_temporalImpulseResponseComplex;
}


double Ganglion::spatialImpulseResponse() const
{
    return m_spatialImpulseResponse;
}

double Ganglion::spatialImpulseResponseComplex() const
{
    return m_spatialImpulseResponseComplex;
}

double Ganglion::temporalImpulseResponse() const
{
    return m_temporalImpulseResponse;
}

double Ganglion::temporalImpulseResponseComplex() const
{
    return m_temporalImpulseResponseComplex;
}







