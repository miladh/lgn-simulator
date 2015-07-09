#ifndef GANGLIONDOG_H
#define GANGLIONDOG_H

#include "ganglion.h"
#include "../math/dog.h"

class GanglionDOG : public Ganglion
{
public:
    GanglionDOG(DOG * dog);
    ~GanglionDOG();

    // Ganglion interface
public:
    void setSpatialImpulseResponseReal(vec2 r);
    void setSpatialImpulseResponseComplex(vec2 k);
    void setTemporalImpulseResponseReal(double t);
    void setTemporalImpulseResponseComplex(double w);


private:
    DOG * m_dog;
};

#endif // GANGLIONDOG_H
