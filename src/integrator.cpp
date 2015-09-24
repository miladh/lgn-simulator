#include "integrator.h"

Integrator::Integrator(IntegratorSettings *settings)
    : m_settings(settings)
    , m_nPointsTemporal(settings->nPointsTemporal())
    , m_nPointsSpatial(settings->nPointsSpatial())
    , m_maxT(settings->maxTime())
    , m_dt(m_maxT/m_nPointsTemporal)
    , m_dw(1./m_nPointsTemporal)
    , m_temporalSamplingFreq(m_nPointsTemporal/m_maxT)
    , m_ds(1.0/m_nPointsSpatial)
    , m_dk(1./m_nPointsSpatial)
    , m_spatialSamplingFreq(m_nPointsSpatial)
{
    //Temporal Grid
    double Nt_2 = ceil(m_nPointsTemporal/2);
    m_timeVec = linspace(0, m_maxT-m_dt, m_nPointsTemporal);
    vec w1 = linspace(0, Nt_2-1, Nt_2);
    vec w2 = linspace(-Nt_2, -w1[1], Nt_2);
    m_temporalFreqs = join_cols(w1,w2)* 1./m_maxT * 2 *PI;

    //Spatial Grid
    double Ns_2 = ceil(m_nPointsSpatial/2);
    m_coordinateVec = linspace(0., 1.0-m_ds, m_nPointsSpatial);
    vec k1 = linspace(0, Ns_2-1, Ns_2);
    vec k2 = linspace(-Ns_2, -k1[1], Ns_2);
    m_spatialFreqs = join_cols(k1,k2) * 2 *PI;

    //initialization
}

Integrator::~Integrator()
{

}

cx_cube Integrator::integrate(cx_cube data)
{
    cx_cube fftData = 0 * data;
    int size[3] = {int(data.n_slices), int(data.n_cols), int(data.n_rows)};

    fftw_complex* in = reinterpret_cast<fftw_complex*> (data.memptr());
    fftw_complex* out = reinterpret_cast<fftw_complex*> (fftData.memptr());
    fftw_plan plan = fftw_plan_dft(3, size, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(plan);
    fftw_destroy_plan(plan);

    return fftData;
}
vec Integrator::timeVec() const
{
    return m_timeVec;
}
vec Integrator::coordinateVec() const
{
    return m_coordinateVec;
}

vec Integrator::temporalFreqVec() const
{
    return m_temporalFreqs;
}

vec Integrator::spatialFreqVec() const
{
    return m_spatialFreqs;
}


int Integrator::nPointsTemporal() const
{
    return m_nPointsTemporal;
}
int Integrator::nPointsSpatial() const
{
    return m_nPointsSpatial;
}











