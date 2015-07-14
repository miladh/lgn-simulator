#include "relay.h"

Relay::Relay(const Config *cfg, Ganglion *ganglion, Stimuli * stim)
    : m_ganglion(ganglion)
    , m_stim(stim)
{
    const Setting & root = cfg->getRoot();
    const Setting &gridLimits = root["gridSettings"]["grid"];
    const Setting &integrationDomain = root["gridSettings"]["integrationDomain"];

    vec3 grid;
    for(int i = 0; i < 3; i++){
        grid[i] = gridLimits[i];
        m_domain[i] = integrationDomain[i];
    }

    m_mesh = linspace(grid[0], grid[1], grid[2]);
    m_response = zeros(grid[2], grid[2]);
    m_responseComplex = zeros(grid[2], grid[2]);
    m_impulseResponse = zeros(grid[2], grid[2]);
    m_impulseResponseComplex = zeros(grid[2], grid[2]);
}

Relay::~Relay()
{
}

void Relay::computeResponse(double t)
{
    mat stim = 0*m_response;
    m_response = 0*m_response;
    m_impulseResponse  = 0* m_impulseResponse;

    double *w = new double [int(m_domain[2])];
    double *x = new double [int(m_domain[2])];
    gauleg(m_domain[0], m_domain[1], x, w, m_domain[2]);

    for(int i = 0; i < int(m_mesh.n_elem); i++){
        for(int j = 0; j < int(m_mesh.n_elem); j++){

            stim(i,j) = m_stim->real({m_mesh[i], m_mesh[j]}, t);
            for(int m = 0; m < int(m_domain[2]); m++){
                for(int n = 0; n < int(m_domain[2]); n++){

                    double T = transferFunctionComplex({x[m], x[n]}, m_stim->w());
                    double Gg = m_ganglion->impulseResponseComplex({x[m], x[n]}, m_stim->w());
                    double s = m_stim->complex({x[m], x[n]}, m_stim->w());

                    double dGr = T * Gg * w[m] * w[n] *
                            cos(m_mesh[i]*x[m]+ m_mesh[j]*x[n] - m_stim->w() * t);
                    m_impulseResponse(i, j) +=  dGr;
                    m_response(i,j) += dGr * s;
                }
            }
        }
    }

    m_stim->setReal(stim);
}

void Relay::computeResponseComplex(double w)
{
    for(int i = 0; i < int(m_mesh.n_elem); i++){
        for(int j = 0; j < int(m_mesh.n_elem); j++){

            double G = transferFunctionComplex({m_mesh[i], m_mesh[j]}, w);
            double s = m_stim->complex({m_mesh[i], m_mesh[j]}, w);

            m_responseComplex(i,j) = G*s;
        }
    }
}


mat Relay::response() const
{
    return m_response;
}

mat Relay::responseComplex() const
{
    return m_responseComplex;
}

mat Relay::impulseRespons() const
{
    return m_impulseResponse;
}

mat Relay::impulseResponsComplex() const
{
    return m_impulseResponseComplex;
}








