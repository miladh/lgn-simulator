#include "response.h"

Response::Response(ImpulseResponse impResFunc, Stimuli stim, Integrator *integrator,
                   vec3 realGrid, vec3 fourierGrid)
    : m_impResFunc(impResFunc)
    , m_stim(stim)
    , m_Integrator(integrator)
    , m_xPoints(linspace(realGrid[0], realGrid[1], realGrid[2]))
    , m_yPoints(linspace(realGrid[0], realGrid[1], realGrid[2]))
    , m_kxPoints(linspace(fourierGrid[0], fourierGrid[1], fourierGrid[2]))
    , m_kyPoints(linspace(fourierGrid[0], fourierGrid[1], fourierGrid[2]))
    , m_response(zeros(realGrid[2], realGrid[2]))
    , m_responseFT(zeros(realGrid[2], realGrid[2]))

{

    computeFT();
}

Response::~Response()
{

}

void Response::computeFT()
{
    for(int i = 0; i < int(m_kxPoints.n_elem); i++){
        for(int j = 0; j < int(m_kyPoints.n_elem); j++){
            double G = m_impResFunc.edogImpulseResponseFunctionFT(
                        m_kxPoints[i], m_kxPoints[j],0.);
            double s = m_stim.patchGratingFT(m_kxPoints[i], m_kxPoints[j], 0);

            m_responseFT(i,j) = G*s;
        }
    }

}

void Response::compute()
{

}

mat Response::realSpace() const
{
    return m_response;
}

mat Response::fourierSpace() const
{
    return m_responseFT;
}

