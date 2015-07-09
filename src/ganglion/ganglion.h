#ifndef GANGLION_H
#define GANGLION_H

#include <armadillo>

using namespace std;
using namespace arma;


class Ganglion
{
public:
    Ganglion();
    ~Ganglion();



    double impulseResponseReal();
    double impulseResponseComplex();

    double spatialImpulseResponseReal() const;
    virtual void setSpatialImpulseResponseReal(vec2 r) = 0;

    double spatialImpulseResponseComplex() const;
    virtual void setSpatialImpulseResponseComplex(vec2 k) = 0;

    double temporalImpulseResponseReal() const;
    virtual void setTemporalImpulseResponseReal(double t) = 0;

    double temporalImpulseResponseComplex() const;
    virtual void setTemporalImpulseResponseComplex(double w) = 0;


protected:
    double m_spatialImpulseResponseReal;
    double m_spatialImpulseResponseComplex;
    double m_temporalImpulseResponseReal;
    double m_temporalImpulseResponseComplex;
};

#endif // GANGLION_H
