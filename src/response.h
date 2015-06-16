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
    Response(ImpulseResponse impResFunc, Stimuli stim,
             Integrator* integrator, vec3 realGrid, vec3 fourierGrid);
    ~Response();
    mat realSpace() const;
    mat fourierSpace() const;

private:
    void computeFT();
    void compute();

    ImpulseResponse m_impResFunc;
    Stimuli m_stim;
    Integrator* m_Integrator;

    vec m_xPoints, m_yPoints;
    vec m_kxPoints, m_kyPoints;
    mat m_response, m_responseFT;
};

#endif // RESPONSE_H
