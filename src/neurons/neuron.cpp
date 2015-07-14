#include "neuron.h"

Neuron::Neuron(const Config *cfg)
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

void Neuron::addGanglionCell(Neuron *neuron,
                                TemporalKernel *tKernel,
                                SpatialKernel *sKernel)
{
    m_ganglionCells.emplace_back(Input{neuron, tKernel, sKernel});
}

void Neuron::addInterNeuron(Neuron *neuron,
                               TemporalKernel *tKernel,
                               SpatialKernel *sKernel)
{
    m_interNeurons.emplace_back(Input{neuron, tKernel, sKernel});
}

void Neuron::addCorticalNeuron(Neuron *neuron,
                                  TemporalKernel *tKernel,
                                  SpatialKernel *sKernel)
{
    m_corticalNeurons.emplace_back(Input{neuron, tKernel, sKernel});
}



vector<Neuron::Input> Neuron::ganglionCells() const
{
    return m_ganglionCells;
}

vector<Neuron::Input> Neuron::relayCells() const
{
    return m_corticalNeurons;
}


vector<Neuron::Input> Neuron::interNeurons() const
{
    return m_interNeurons;
}

vector<Neuron::Input> Neuron::corticalNeurons() const
{
    return m_corticalNeurons;
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










