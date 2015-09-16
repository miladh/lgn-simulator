#include "stimuli.h"


Stimuli::Stimuli(const Config *cfg)
{
    const Setting & root = cfg->getRoot();
    m_w = root["stimuliSettings"]["w"];
    m_k[0] = root["stimuliSettings"]["kx"];
    m_k[1] = root["stimuliSettings"]["ky"];
    m_nPoints = root["spatialDomainSettings"]["nPoints"];
    m_nSteps = root["dynamicSettings"]["nSteps"];

    m_spatial = zeros<cube>(m_nPoints, m_nPoints, m_nSteps);
    m_frequency = zeros<cx_cube>(m_nPoints, m_nPoints, m_nSteps);


    //Spatial Mesh
    m_spatialMesh = linspace(-0.5, 0.5, m_nPoints);
    double dr = m_spatialMesh(1) - m_spatialMesh(0);
    double N_2 = ceil(m_nPoints/2.);
    double df = 1./dr/m_nPoints;
    double fs = 1./dr;

    m_spatialFreqs = linspace(-N_2*df, (m_nPoints - 1. - N_2)*df, m_nPoints);
    m_spatialFreqs*=2*PI;

    //Temporal Mesh
    m_temporalMesh = linspace(-0.5, 0.5, m_nSteps);
    double dt = m_temporalMesh(1) - m_temporalMesh(0);
    double Nt_2 = ceil(m_nSteps/2.);
    double df_t = 1./dt/m_nSteps;
    double fs_t = 1./dt;

    m_temporalFreqs = linspace(-Nt_2*df_t, (m_nSteps - 1. - Nt_2)*df_t, m_nSteps);
    m_temporalFreqs*=2*PI;
}

Stimuli::~Stimuli()
{

}

void Stimuli::computeSpatial()
{
    for(int k = 0; k < m_nSteps; k++){
        for(int i = 0; i < m_nPoints; i++){
            for(int j = 0; j < m_nPoints; j++){
                m_spatial(i,j,k) = spatial({m_spatialMesh[i], m_spatialMesh[j]},
                                           m_temporalMesh[k]);
            }
        }
    }

}


void Stimuli::computeFrequency()
{
    for(int k = 0; k < m_nSteps; k++){
        for(int i = 0; i < m_nPoints; i++){
            for(int j = 0; j < m_nPoints; j++){
                m_frequency(i,j,k) = frequency({m_spatialFreqs[i],
                                              m_spatialFreqs[j]},
                                             m_temporalFreqs[k]);

            }
        }
    }
}


double Stimuli::w() const
{
    return m_w;
}

void Stimuli::setSpatial(cube spatial)
{
    m_spatial = spatial;
}


cube Stimuli::spatial() const
{
    return m_spatial;
}

cx_cube Stimuli::frequency() const
{
    return m_frequency;
}

