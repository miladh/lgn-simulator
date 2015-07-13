#ifndef GANGLIONDOG_H
#define GANGLIONDOG_H

#include "ganglion.h"
#include "../math/dog.h"

class GanglionDOG : public Ganglion
{
public:
    GanglionDOG(const Config *cfg);
    ~GanglionDOG();

    // Ganglion interface
public:
    void setSpatialImpulseResponse(vec2 r);
    void setSpatialImpulseResponseComplex(vec2 k);
    void setTemporalImpulseResponse(double t);
    void setTemporalImpulseResponseComplex(double w);


private:
    DOG * m_dog;
    double m_tau_rg = 0.0;  //[s]
};

#endif // GANGLIONDOG_H
