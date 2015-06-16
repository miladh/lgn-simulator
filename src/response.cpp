#include "response.h"

Response::Response(ImpulseResponse impResFunc, Stimuli stim, Integrator *integrator,
                   vec3 grid)
    : m_impResFunc(impResFunc)
    , m_stim(stim)
    , m_Integrator(integrator)
    , m_xPoints(linspace(grid[0], grid[1], grid[2]))
    , m_yPoints(linspace(grid[0], grid[1], grid[2]))
    , m_response(zeros(grid[2], grid[2]))
{

    compute();
}

Response::~Response()
{

}

void Response::compute()
{
    for(int  i= 0; i < int(m_xPoints.n_elem); i++){
        for(int  j= 0; j < int(m_yPoints.n_elem); j++){
            m_response(i,j) = i+j;
        }


    }

}
mat Response::response() const
{
    return m_response;
}

