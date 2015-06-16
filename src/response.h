#ifndef RESPONSE_H
#define RESPONSE_H

#include <armadillo>
#include <impulseResponse.h>
#include <stimuli.h>
#include <integrator.h>

using namespace arma;
using namespace std;

class Response
{
public:
    Response(ImpulseResponse impResFunc, Stimuli stim, Integrator* integrator,
             vec3 grid);
    ~Response();
    mat response() const;

private:
    void compute();

    ImpulseResponse m_impResFunc;
    Stimuli m_stim;
    Integrator* m_Integrator;

    vec m_xPoints;
    vec m_yPoints;
    mat m_response;
};

#endif // RESPONSE_H
