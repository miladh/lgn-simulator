#include "integrator.h"

using namespace lgnSimulator;


Integrator::Integrator(const int nt,
                       const double temporalResolution,
                       const int ns,
                       const double spatialResolution)

    : m_nPointsTemporal(pow(2,nt))
    , m_nPointsSpatial(pow(2,ns))
    , m_temporalResolution(temporalResolution)
    , m_spatialResolution(spatialResolution)
    , m_timeInterval(m_nPointsTemporal  * m_temporalResolution)
    , m_lengthInterval(m_nPointsSpatial * m_spatialResolution)
    , m_temporalSamplingFreq(2.* core::pi / m_temporalResolution)
    , m_spatialSamplingFreq(2.* core::pi / m_spatialResolution)
    , m_temporalFreqResolution(m_temporalSamplingFreq / m_nPointsTemporal)
    , m_spatialFreqResolution(m_spatialSamplingFreq / m_nPointsSpatial)
{
    //Temporal Grid
    m_temporalFreqs = FFTHelper::fftFreq(m_nPointsTemporal, m_temporalResolution)*2*core::pi;
    m_timeVec = linspace(0,m_nPointsTemporal-1,m_nPointsTemporal)
                * m_temporalResolution;

    //Spatial Grid
    m_spatialFreqs = FFTHelper::fftFreq(m_nPointsSpatial, m_spatialResolution)*2*core::pi;
    m_spatialVec = linspace(-m_nPointsSpatial/2,(m_nPointsSpatial-1)/2, m_nPointsSpatial)
                  * m_spatialResolution;

}

Integrator::~Integrator()
{
}

cx_cube Integrator::backwardFFT(cx_cube data) const
{
    cx_cube fftData = 0 * data;
    int size[3] = {int(data.n_slices), int(data.n_cols), int(data.n_rows)};

    fftw_complex* in = reinterpret_cast<fftw_complex*> (data.memptr());
    fftw_complex* out = reinterpret_cast<fftw_complex*> (fftData.memptr());
    fftw_plan plan = fftw_plan_dft(3, size, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(plan);
    fftw_destroy_plan(plan);

    fftData *= m_temporalFreqResolution
            * m_spatialFreqResolution
            * m_spatialFreqResolution
            /8./core::pi/core::pi/core::pi;


    //fftShift
    for(int i = 0; i < int(fftData.n_slices); i++){
        fftData.slice(i) = FFTHelper::fftShift(fftData.slice(i));
    }

    return fftData;
}



cx_cube Integrator::forwardFFT(cx_cube data) const
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


cx_mat Integrator::backwardFFT(cx_mat data) const
{
    cx_mat fftData = 0 * data;
    int size[2] = {int(data.n_cols), int(data.n_rows)};

    fftw_complex* in = reinterpret_cast<fftw_complex*> (data.memptr());
    fftw_complex* out = reinterpret_cast<fftw_complex*> (fftData.memptr());
    fftw_plan plan = fftw_plan_dft(2, size, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(plan);
    fftw_destroy_plan(plan);

    fftData *= m_spatialFreqResolution*m_spatialFreqResolution /4./core::pi/core::pi;

    //fftShift
    fftData = FFTHelper::fftShift(fftData);

    return fftData;
}


cx_mat Integrator::forwardFFT(cx_mat data) const
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
    return m_temporalResolution;
}
double Integrator::spatialResolution() const
{
    return m_spatialResolution;
}
double Integrator::temporalFreqResolution() const
{
    return m_temporalFreqResolution;
}
double Integrator::spatialFreqResolution() const
{
    return m_spatialFreqResolution;
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


Integrator createIntegrator(const YAML::Node &cfg)
{
    int nt = cfg["nt"].as<int>();
    int ns = cfg["ns"].as<int>();
    double dt = cfg["dt"].as<double>();
    double ds = cfg["ds"].as<double>();

    return Integrator(nt, dt, ns, ds);

}
