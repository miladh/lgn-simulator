#include "integrator.h"

Integrator::Integrator(IntegratorSettings *settings)
    : m_settings(settings)
    , m_nPointsTemporal(settings->nPointsTemporal())
    , m_nPointsSpatial(settings->nPointsSpatial())
    , m_dt(settings->temporalResolution())
    , m_ds(settings->spatialResolution())
    , m_dw(2.* PI/m_nPointsTemporal/m_dt)
    , m_dk(2.* PI/m_nPointsSpatial/m_ds)
    //    , m_temporalSamplingFreq(m_nPointsTemporal/m_maxT)
    //    , m_spatialSamplingFreq(m_nPointsSpatial)
{
    //Temporal Grid
//    m_timeVec = linspace(0, m_nPointsTemporal-1 , m_nPointsTemporal)*m_dt;
    m_timeVec = linspace(-m_nPointsTemporal/2,
                               m_nPointsTemporal/2-1,
                               m_nPointsTemporal)*m_dt;

    //    m_temporalFreqs = FFTHelper::fftFreq(m_nPointsTemporal, m_dt)*2*PI;
    m_temporalFreqs = linspace(-m_nPointsTemporal/2,
                               m_nPointsTemporal/2-1,
                               m_nPointsTemporal)*m_dw;

    //Spatial Grid
    m_coordinateVec = linspace(-m_nPointsSpatial/2,
                               m_nPointsSpatial/2-1,
                               m_nPointsSpatial)*m_ds;
    //    m_spatialFreqs = FFTHelper::fftFreq(m_nPointsSpatial, m_ds)*2*PI;
    m_spatialFreqs = linspace(-m_nPointsSpatial/2,
                              m_nPointsSpatial/2-1,
                              m_nPointsSpatial)*m_dk;

    cout << "max x:" << m_coordinateVec(m_nPointsSpatial-1) << endl;
    cout << "dw :" << m_dw << endl;
    cout << "dk:" << m_dk << endl;
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

    for(int k=0; k < int(data.n_slices); k++){
        for(int i=0; i < int(data.n_rows); i++){
            for(int j=0; j < int(data.n_cols); j++){
                fftData.slice(k)(i,j) *= pow(-1,i+j-k);
            }
        }
    }


    fftData *= m_dw * m_dk * m_dk /8./PI/PI/PI;
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


    for(int i=0; i < int(data.n_rows); i++){
        for(int j=0; j < int(data.n_cols); j++){
            fftData(i,j) *= pow(-1,i+j);
        }
    }

    fftData *= m_dk * m_dk /4./PI/PI;
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













