#ifndef GANGLION_H
#define GANGLION_H

#include <libconfig.h++>
#include <armadillo>


using namespace std;
using namespace libconfig;
using namespace arma;


class Ganglion
{
public:
    Ganglion(const Config *cfg);
    ~Ganglion();



    double impulseResponse(vec2 r, double t);
    double impulseResponseComplex(vec2 k, double w);

    double spatialImpulseResponse() const;
    virtual void setSpatialImpulseResponse(vec2 r) = 0;

    double spatialImpulseResponseComplex() const;
    virtual void setSpatialImpulseResponseComplex(vec2 k) = 0;

    double temporalImpulseResponse() const;
    virtual void setTemporalImpulseResponse(double t) = 0;

    double temporalImpulseResponseComplex() const;
    virtual void setTemporalImpulseResponseComplex(double w) = 0;


protected:
    double m_spatialImpulseResponse = 0;
    double m_spatialImpulseResponseComplex = 0;
    double m_temporalImpulseResponse = 0;
    double m_temporalImpulseResponseComplex = 0;
};

#endif // GANGLION_H
