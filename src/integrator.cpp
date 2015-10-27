#include "integrator.h"

Integrator::Integrator(int nt, double dt, int ns, double ds)
    : m_nPointsTemporal(pow(2,nt))
    , m_dt(dt)
    , m_nPointsSpatial(pow(2,ns))
    , m_ds(ds)
    , m_dw(2.* PI/m_nPointsTemporal/m_dt)
    , m_dk(2.* PI/m_nPointsSpatial/m_ds)
{
    //Temporal Grid
    m_timeVec = linspace(0, m_nPointsTemporal-1 , m_nPointsTemporal)*m_dt;
    m_temporalFreqs = FFTHelper::fftFreq(m_nPointsTemporal, m_dt)*2*PI;


    //Spatial Grid
    m_coordinateVec = linspace(-m_nPointsSpatial/2, m_nPointsSpatial/2-1,
                               m_nPointsSpatial)*m_ds;
    m_spatialFreqs = FFTHelper::fftFreq(m_nPointsSpatial, m_ds)*2*PI;


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

    fftData *= m_dw * m_dk * m_dk /8./PI/PI/PI;


    //fftShift
    for(int i = 0; i < int(fftData.n_slices); i++){
        fftData.slice(i) = FFTHelper::fftShift(fftData.slice(i));
    }

    return fftData;
}

cx_mat Integrator::integrate(cx_mat data)
{
    cx_mat fftData = 0 * data;
    int size[2] = {int(data.n_cols), int(data.n_rows)};

    fftw_complex* in = reinterpret_cast<fftw_complex*> (data.memptr());
    fftw_complex* out = reinterpret_cast<fftw_complex*> (fftData.memptr());
    fftw_plan plan = fftw_plan_dft(2, size, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(plan);
    fftw_destroy_plan(plan);

    fftData *= m_dk * m_dk /4./PI/PI;

    //fftShift
    fftData = FFTHelper::fftShift(fftData);

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
double Integrator::dt() const
{
    return m_dt;
}
double Integrator::ds() const
{
    return m_ds;
}
double Integrator::temporalFreqResolution() const
{
    return m_dw;
}
double Integrator::spatialFreqResolution() const
{
    return m_dk;
}



Integrator createIntegrator(const Config *cfg)
{
    const Setting & root = cfg->getRoot();
    int ns = root["integratorSettings"]["ns"];
    int nt = root["integratorSettings"]["nt"];
    double dt = root["integratorSettings"]["dt"];
    double ds = root["integratorSettings"]["ds"];

    return Integrator(nt, dt, ns, ds);

}
