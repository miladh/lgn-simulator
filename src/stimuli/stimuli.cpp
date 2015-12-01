#include "stimuli.h"


Stimulus::Stimulus(Integrator *integrator)
    : m_integrator(integrator)
{
    int nPointsTemporal = integrator->nPointsTemporal();
    int nPointsSpatial = integrator->nPointsSpatial();

    m_spatioTemporal = zeros<cube>(nPointsSpatial, nPointsSpatial, nPointsTemporal);
    m_fourierTransform = zeros<cx_cube>(nPointsSpatial, nPointsSpatial, nPointsTemporal);

    //Temporal Mesh
    m_timeVec = integrator->timeVec();
    m_temporalFreqs = integrator->temporalFreqVec();

    //Spatial Mesh
    m_coordinateVec = integrator->coordinateVec();
    m_spatialFreqs = integrator->spatialFreqVec();


}

Stimulus::~Stimulus()
{
}

cube Stimulus::spatioTemporal() const
{
    return m_spatioTemporal;
}

cx_cube Stimulus::fourierTransform() const
{
    return m_fourierTransform;
}

void Stimulus::computeSpatiotemporalAnalytic()
{
    for(int k = 0; k < int(m_spatioTemporal.n_slices); k++){
        for(int i = 0; i < int(m_spatioTemporal.n_rows); i++){
            for(int j = 0; j < int(m_spatioTemporal.n_cols); j++){
                m_spatioTemporal(i,j,k) = valueAtPoint({m_coordinateVec[i],
                                                        m_coordinateVec[j]},
                                                       m_timeVec[k]);
            }
        }
    }
}

void Stimulus::computeFourierTransformAnalytic()
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


void Stimulus::clearSpatioTemporal()
{
    m_spatioTemporal.clear();
}


void Stimulus::clearFourierTransform()
{

    m_fourierTransform.clear();

}
