#ifndef RESPONSE_H
#define RESPONSE_H

#include <armadillo>
#include <impulseResponse.h>
#include <stimuli.h>
#include <integrator.h>

#include <lib.h>

using namespace arma;
using namespace std;

class Response
{
public:
    Response(ImpulseResponse impResFunc, Stimuli stim,
             Integrator* integrator, vec3 realGrid, vec3 complexGrid, vec3 integrationDomain);
    ~Response();
    mat real() const;
    mat complex() const;

    void compute(double t);
    void computeComplex(double w);

private:

    ImpulseResponse m_impResFunc;
    Stimuli m_stim;
    Integrator* m_Integrator;

    vec m_xPoints, m_yPoints;
    vec m_kxPoints, m_kyPoints;
    vec m_domain;
    mat m_response, m_responseComplex;
};

#endif // RESPONSE_H
