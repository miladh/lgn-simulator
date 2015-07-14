#include "corticalcell.h"

CorticalCell::CorticalCell(const Config * cfg)
    : Neuron(cfg)
{

}

CorticalCell::~CorticalCell()
{

}

void CorticalCell::computeResponse()
{

    m_response = 0 * m_response;
    m_impulseResponse = 0 * m_impulseResponse;
}




