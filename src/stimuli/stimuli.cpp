#include "stimuli.h"


Stimuli::Stimuli(const Config *cfg, Integrator integrator)
    : m_integrator(integrator)
{
    const Setting & root = cfg->getRoot();
    m_k[0] = root["stimuliSettings"]["kx"];
    m_k[1] = root["stimuliSettings"]["ky"];

    int nPointsTemporal = integrator.nPointsTemporal();
    int nPointsSpatial = integrator.nPointsSpatial();

    m_spatioTemporal = zeros<cube>(nPointsSpatial, nPointsSpatial, nPointsTemporal);
    m_fourierTransform = zeros<cx_cube>(nPointsSpatial, nPointsSpatial, nPointsTemporal);

    //Temporal Mesh
    timeVec = integrator.timeVec();
    m_temporalFreqs = integrator.temporalFreqVec();

    //Spatial Mesh
    m_coordinateVec = integrator.coordinateVec();
    m_spatialFreqs =integrator.spatialFreqVec();


    m_w = m_temporalFreqs[m_temporalFreqs.n_elem/2];
    m_w = m_temporalFreqs[0];
    m_k[0] = m_spatialFreqs[5];
    cout <<"m_w: "<< m_w << endl;
    cout <<"m_kx: "<< m_k[0] << endl;

}

Stimuli::~Stimuli()
{

}

void Stimuli::computeSpatiotemporal()
{
    for(int k = 0; k < int(m_spatioTemporal.n_slices); k++){
        for(int i = 0; i < int(m_spatioTemporal.n_rows); i++){
            for(int j = 0; j < int(m_spatioTemporal.n_cols); j++){
                m_spatioTemporal(i,j,k) = valueAtPoint({m_coordinateVec[i],
                                                        m_coordinateVec[j]},
                                                       timeVec[k]);
            }
        }
    }

}


void Stimuli::computeFourierTransform()
{

    for(int k = 0; k < int(m_fourierTransform.n_slices); k++){
        for(int i = 0; i < int(m_fourierTransform.n_rows); i++){
            for(int j = 0; j < int(m_fourierTransform.n_cols); j++){
                m_fourierTransform(i,j,k) =
                        fourierTransformAtFrequency({m_spatialFreqs[i],
                                                     m_spatialFreqs[j]},
                                                    m_temporalFreqs[k]);

            }
        }
    }
}
cube Stimuli::spatioTemporal() const
{
    return m_spatioTemporal;
}

cx_cube Stimuli::fourierTransform() const
{
    return m_fourierTransform;
}



