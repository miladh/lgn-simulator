#include "response.h"

Response::Response(ImpulseResponse impResFunc, Stimuli stim, Integrator *integrator,
                   vec3 realGrid, vec3 complexGrid, vec3 integrationDomain)
    : m_impResFunc(impResFunc)
    , m_stim(stim)
    , m_Integrator(integrator)
    , m_xPoints(linspace(realGrid[0], realGrid[1], realGrid[2]))
    , m_yPoints(linspace(realGrid[0], realGrid[1], realGrid[2]))
    , m_kxPoints(linspace(complexGrid[0], complexGrid[1], complexGrid[2]))
    , m_kyPoints(linspace(complexGrid[0], complexGrid[1], complexGrid[2]))
    , m_domain(integrationDomain)
    , m_response(ones(realGrid[2], realGrid[2]))
    , m_responseComplex(zeros(complexGrid[2], complexGrid[2]))

{

}

Response::~Response()
{

}

void Response::computeComplex(double w)
{

    for(int i = 0; i < int(m_kxPoints.n_elem); i++){
        for(int j = 0; j < int(m_kyPoints.n_elem); j++){
            double G = m_impResFunc.edogComplex(
                        m_kxPoints[i], m_kxPoints[j], w);
            double s = m_stim.patchGratingComplex(m_kxPoints[i], m_kxPoints[j], w);

            m_responseComplex(i,j) = G*s;
        }
    }

}

void Response::compute(double t)
{

    m_response = 0*m_response;
    double *w = new double [int(m_domain[2])];
    double *x = new double [int(m_domain[2])];
    gauleg(m_domain[0], m_domain[1], x,w, m_domain[2]);


    for(int i = 0; i < int(m_xPoints.n_elem); i++){
        for(int j = 0; j < int(m_yPoints.n_elem); j++){

            for(int m = 0; m < int(m_domain[2]); m++){
                for(int n = 0; n < int(m_domain[2]); n++){
                    for(int o = 0; o < int(m_domain[2]); o++){

                        double G = m_impResFunc.edogComplex(
                                    x[m], x[n], x[o]);
                        double s = m_stim.patchGratingComplex(x[m], x[n], x[o]);

                        m_response(i,j) += s * w[m] * w[n] * w[o]
                                * cos(m_xPoints[i]*x[m]+ m_yPoints[j]*x[n] - x[o] * t);
                    }
                }
            }
        }
    }

}

mat Response::real() const
{
    return m_response;
}

mat Response::complex() const
{
    return m_responseComplex;
}

