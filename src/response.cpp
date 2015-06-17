#include "response.h"

Response::Response(ImpulseResponse impResFunc, Stimuli stim, Integrator *integrator,
                   vec3 realGrid, vec3 complexGrid)
    : m_impResFunc(impResFunc)
    , m_stim(stim)
    , m_Integrator(integrator)
    , m_xPoints(linspace(realGrid[0], realGrid[1], realGrid[2]))
    , m_yPoints(linspace(realGrid[0], realGrid[1], realGrid[2]))
    , m_kxPoints(linspace(complexGrid[0], complexGrid[1], complexGrid[2]))
    , m_kyPoints(linspace(complexGrid[0], complexGrid[1], complexGrid[2]))
    , m_response(ones(realGrid[2], realGrid[2]))
    , m_responseFT(zeros(realGrid[2], realGrid[2]))

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

            m_responseFT(i,j) = G*s;
        }
    }

}

void Response::compute()
{

}

mat Response::real() const
{
    return m_response;
}

mat Response::complex() const
{
    return m_responseFT;
}

