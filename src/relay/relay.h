#ifndef RELAY_H
#define RELAY_H

#include <libconfig.h++>
#include <armadillo>
#include "../ganglion/ganglion.h"
#include "../ganglion/gangliondog.h"
#include "../stimuli/stimuli.h"
#include "../stimuli/patchgrating.h"

#include <lib.h>

using namespace std;
using namespace libconfig;
using namespace arma;


class Relay
{
public:
    Relay(const Config * cfg, Ganglion * ganglion, Stimuli * stim);
    ~Relay();

    virtual double transferFunctionComplex(vec2 kVec, double w) = 0;
    void computeResponse(double t);
    void computeResponseComplex(double w);

    mat response() const;
    mat responseComplex() const;
    mat impulseRespons() const;
    mat impulseResponsComplex() const;



protected:
    mat m_response, m_responseComplex;
    mat m_impulseResponse, m_impulseResponseComplex;


    Ganglion * m_ganglion;
    Stimuli * m_stim;

    vec m_mesh;
    vec3 m_domain;

};

#endif // RELAY_H
