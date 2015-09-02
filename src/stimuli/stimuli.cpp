#include "stimuli.h"


Stimuli::Stimuli(const Config *cfg)
{
    const Setting & root = cfg->getRoot();
    m_w = root["stimuliSettings"]["w"];
    m_k[0] = root["stimuliSettings"]["kx"];
    m_k[1] = root["stimuliSettings"]["ky"];
    m_nPoints = root["spatialDomainSettings"]["nPoints"];
    m_spatialMesh = linspace(-0.5, 0.5, m_nPoints);


    m_spatial = zeros<mat>(m_nPoints, m_nPoints);
    m_frequency = zeros<cx_mat>(m_nPoints, m_nPoints);

}

Stimuli::~Stimuli()
{

}

void Stimuli::computeSpatial(double t)
{

    for(int i = 0; i < m_nPoints; i++){
        for(int j = 0; j < m_nPoints; j++){
            m_spatial(i,j) = spatial({m_spatialMesh[i], m_spatialMesh[j]}, t);
        }
    }
}

void Stimuli::computeFrequency()
{
    m_frequency = fft2(spatial());
}


double Stimuli::w() const
{
    return m_w;
}

void Stimuli::setSpatial(mat spatial)
{
    m_spatial = spatial;
}


mat Stimuli::spatial() const
{
    return m_spatial;
}

cx_mat Stimuli::frequency() const
{
    return m_frequency;
}

