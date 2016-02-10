#include "integrator.h"

using namespace lgnSimulator;


Integrator::Integrator(int nt, double dt, int ns, double ds)
    : m_nPointsTemporal(pow(2,nt))
    , m_nPointsSpatial(pow(2,ns))
    , m_dt(dt)
    , m_ds(ds)
    , m_timeInterval(m_nPointsTemporal  * m_dt)
    , m_lengthInterval(m_nPointsSpatial * m_ds)
    , m_temporalSamplingFreq(2.* PI / m_dt)
    , m_spatialSamplingFreq(2.* PI / m_ds)
    , m_dw(m_temporalSamplingFreq / m_nPointsSpatial)
    , m_dk(m_spatialSamplingFreq / m_nPointsSpatial)
{
    //Temporal Grid
    m_timeVec = linspace(0, m_nPointsTemporal-1 , m_nPointsTemporal)*m_dt;
    m_temporalFreqs = FFTHelper::fftFreq(m_nPointsTemporal, m_dt)*2*PI;

    //Spatial Grid
    m_spatialVec = linspace(-m_nPointsSpatial/2 , (m_nPointsSpatial-1)/2, m_nPointsSpatial)*m_ds;
    m_spatialFreqs = FFTHelper::fftFreq(m_nPointsSpatial, m_ds)*2*PI;

}

Integrator::~Integrator()
{
}

cx_cube Integrator::backwardFFT(cx_cube data)
{
    cx_cube fftData = 0 * data;
    int size[3] = {int(data.n_slices), int(data.n_cols), int(data.n_rows)};

    for(int k = 0; k < int(fftData.n_slices); k++){
        for(int i = 0; i < int(fftData.n_cols); i++){
            for(int j = 0; j < int(fftData.n_rows); j++){
                data(i,j,k) *=pow(-1, i+j);
            }
        }
    }


    fftw_complex* in = reinterpret_cast<fftw_complex*> (data.memptr());
    fftw_complex* out = reinterpret_cast<fftw_complex*> (fftData.memptr());
    fftw_plan plan = fftw_plan_dft(3, size, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(plan);
    fftw_destroy_plan(plan);

    fftData *= m_dw * m_dk * m_dk /8./PI/PI/PI;

//    for(int k = 0; k < int(fftData.n_slices); k++){
//        for(int i = 0; i < int(fftData.n_cols); i++){
//            for(int j = 0; j < int(fftData.n_rows); j++){
//                fftData(i,j,k) *=pow(-1, i+j);
//            }
//        }
//    }


//        //fftShift
//        for(int i = 0; i < int(fftData.n_slices); i++){
//            fftData.slice(i) = FFTHelper::fftShift(fftData.slice(i));
//        }


    return fftData;
}



cx_cube Integrator::forwardFFT(cx_cube data)
{
    //fftShift
    for(int i = 0; i < int(data.n_slices); i++){
        data.slice(i) = FFTHelper::fftShift(data.slice(i));
    }

    cx_cube ifftData = 0 * data;
    int size[3] = {int(data.n_slices), int(data.n_cols), int(data.n_rows)};

    fftw_complex* in = reinterpret_cast<fftw_complex*> (data.memptr());
    fftw_complex* out = reinterpret_cast<fftw_complex*> (ifftData.memptr());
    fftw_plan plan = fftw_plan_dft(3, size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_execute(plan);
    fftw_destroy_plan(plan);
    return ifftData;

}


cx_mat Integrator::backwardFFT(cx_mat data)
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


cx_mat Integrator::forwardFFT(cx_mat data)
{
    //fftShift
    data = FFTHelper::fftShift(data);

    cx_mat ifftData = 0 * data;
    int size[2] = {int(data.n_cols), int(data.n_rows)};

    fftw_complex* in = reinterpret_cast<fftw_complex*> (data.memptr());
    fftw_complex* out = reinterpret_cast<fftw_complex*> (ifftData.memptr());
    fftw_plan plan = fftw_plan_dft(2, size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_execute(plan);
    fftw_destroy_plan(plan);

    return ifftData;
}

vec Integrator::timeVec() const
{
    return m_timeVec;
}
vec Integrator::spatialVec() const
{
    return m_spatialVec;
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
double Integrator::temporalResolution() const
{
    return m_dt;
}
double Integrator::spatialResolution() const
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

double Integrator::timeInterval() const
{
    return m_timeInterval;
}

double Integrator::lengthInterval() const
{
    return m_lengthInterval;
}

double Integrator::temporalSamplingFreq() const
{
    return m_temporalSamplingFreq;
}

double Integrator::spatialSamplingFreq() const
{
    return m_spatialSamplingFreq;
}


Integrator createIntegrator(const YAML::Node *cfg)
{
    int nt = (*cfg)["nt"].as<int>();
    int ns =(*cfg)["ns"].as<int>();
    double dt = (*cfg)["dt"].as<double>();
    double ds = (*cfg)["ds"].as<double>();

    return Integrator(nt, dt, ns, ds);

}
