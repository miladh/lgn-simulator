#include "stimuli.h"

using namespace lgnSimulator;


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
    m_spatialVec = integrator->spatialVec();
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


void Stimulus::clearSpatioTemporal()
{
    m_spatioTemporal.clear();
}


void Stimulus::clearFourierTransform()
{

    m_fourierTransform.clear();

}

string Stimulus::type() const
{
    return m_type;
}
