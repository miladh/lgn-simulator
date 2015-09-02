#include "stimuli.h"


Stimuli::Stimuli(const Config *cfg)
{
    const Setting & root = cfg->getRoot();
    m_w = root["stimuliSettings"]["w"];
    m_k[0] = root["stimuliSettings"]["kx"];
    m_k[1] = root["stimuliSettings"]["ky"];


    int nPoints = root["spatialDomainSettings"]["nPoints"];
    double dr = root["spatialDomainSettings"]["dr"];

    m_spatialMesh = linspace(-0.5, 0.5, nPoints);


    m_spatial = zeros<mat>(nPoints, nPoints);
    m_frequency = zeros<cx_mat>(nPoints, nPoints);

}

Stimuli::~Stimuli()
{

}

void Stimuli::computeSpatial(double t)
{

    for(int i = 0; i < int(m_spatialMesh.n_elem); i++){
        for(int j = 0; j < int(m_spatialMesh.n_elem); j++){
            m_spatial(i,j) = spatial({m_spatialMesh[i], m_spatialMesh[j]}, t);
        }
    }

    computeFrequency();

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

