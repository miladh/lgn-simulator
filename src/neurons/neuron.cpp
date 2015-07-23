#include "neuron.h"

Neuron::Neuron(const Config *cfg, Stimuli *stim)
    : m_stim(stim)
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

Neuron::~Neuron()
{

}

void Neuron::computeResponse(double t)
{

    mat stim = 0*m_response;
    m_response = 0*m_response;
    double *w = new double [int(m_domain[2])];
    double *x = new double [int(m_domain[2])];
    gauleg(m_domain[0], m_domain[1], x, w, m_domain[2]);

    for(int i = 0; i < int(m_mesh.n_elem); i++){
        for(int j = 0; j < int(m_mesh.n_elem); j++){

            stim(i,j) = m_stim->real({m_mesh[i], m_mesh[j]}, t);
            for(int m = 0; m < int(m_domain[2]); m++){
                for(int n = 0; n < int(m_domain[2]); n++){

                    double Scomplex = m_stim->complex({x[m], x[n]}, m_stim->w());
                    double Gcomplex = impulseResponseComplex({x[m], x[n]}, m_stim->w());

                    double dGr = Gcomplex * w[m] * w[n] *
                            cos(m_mesh[i]*x[m]+ m_mesh[j]*x[n] - m_stim->w() * t);
                    m_response(i,j) += dGr * Scomplex;

                }
            }
        }
    }
    m_stim->setReal(stim);

}

void Neuron::computeResponseComplex(double w)
{
    cout << "computeResponseComplex: Not implemented!" << endl;
    for(int i = 0; i < int(m_mesh.n_elem); i++){
        for(int j = 0; j < int(m_mesh.n_elem); j++){

            double Gcomplex = impulseResponseComplex({m_mesh[i], m_mesh[j]}, w);
            double Scomplex = m_stim->complex({m_mesh[i], m_mesh[j]}, w);

            m_impulseResponseComplex(i,j) = Gcomplex;
            m_responseComplex(i,j) = Gcomplex*Scomplex;
        }
    }

}

void Neuron::computeImpulseResponse(double t)
{
    m_impulseResponse  = 0* m_impulseResponse;
    double *w = new double [int(m_domain[2])];
    double *x = new double [int(m_domain[2])];
    gauleg(m_domain[0], m_domain[1], x, w, m_domain[2]);

    for(int i = 0; i < int(m_mesh.n_elem); i++){
        for(int j = 0; j < int(m_mesh.n_elem); j++){

            for(int m = 0; m < int(m_domain[2]); m++){
                for(int n = 0; n < int(m_domain[2]); n++){
                    for(int o = 0; o < int(m_domain[2]); o++){
                        m_impulseResponse(i, j) +=
                                impulseResponseComplex({x[m], x[n]},x[o])
                                * w[m] * w[n] * w[o]
                                * cos(m_mesh[i]*x[m]+ m_mesh[j]*x[n] - x[o]* t);
                    }

                }
            }
        }
    }
}

void Neuron::computeImpulseResponseComplex(double w)
{
    cout << "computeImpulseResponseComplex: Not implemented!" << endl;
}



void Neuron::addGanglionCell(Neuron *neuron,
                             SpatialKernel *sKernel,
                             TemporalKernel *tKernel)
{
    m_ganglionCells.emplace_back(Input{neuron, sKernel, tKernel});
}

void Neuron::addRelayCell(Neuron *neuron,
                          SpatialKernel *sKernel,
                          TemporalKernel *tKernel)
{
    m_relayCells.emplace_back(Input{neuron, sKernel, tKernel});
}

void Neuron::addInterNeuron(Neuron *neuron,
                            SpatialKernel *sKernel,
                            TemporalKernel *tKernel)
{
    m_interNeurons.emplace_back(Input{neuron, sKernel, tKernel});
}

void Neuron::addCorticalNeuron(Neuron *neuron,
                               SpatialKernel *sKernel,
                               TemporalKernel *tKernel)
{
    m_corticalNeurons.emplace_back(Input{neuron, sKernel, tKernel});
}


vector<Neuron::Input> Neuron::ganglionCells() const
{
    return m_ganglionCells;
}

vector<Neuron::Input> Neuron::relayCells() const
{
    return m_relayCells;
}


vector<Neuron::Input> Neuron::interNeurons() const
{
    return m_interNeurons;
}

vector<Neuron::Input> Neuron::corticalNeurons() const
{
    return m_corticalNeurons;
}


string Neuron::cellType() const
{
    return m_cellType;
}

mat Neuron::response() const
{
    return m_response;
}

mat Neuron::impulseResponse() const
{
    return m_impulseResponse;
}

mat Neuron::responseComplex() const
{
    return m_responseComplex;
}

mat Neuron::impulseResponseComplex() const
{
    return m_impulseResponseComplex;
}










