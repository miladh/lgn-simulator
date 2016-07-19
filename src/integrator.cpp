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
    m_temporalFreqs = FFTHelper::fftFreq(m_nPointsTemporal,
                                         m_temporalResolution)*-2*core::pi;
    m_timeVec = linspace(0,m_nPointsTemporal-1, m_nPointsTemporal)
            * m_temporalResolution;

    //Spatial Grid
    m_spatialFreqs = FFTHelper::fftFreq(m_nPointsSpatial, m_spatialResolution)*2*core::pi;
    m_spatialVec = linspace(-m_nPointsSpatial/2, m_nPointsSpatial/2-1, m_nPointsSpatial)
            * m_spatialResolution;


    m_real.set_size(m_nPointsSpatial, m_nPointsSpatial, m_nPointsTemporal);
    m_complex.set_size(m_nPointsSpatial, m_nPointsSpatial, m_nPointsTemporal);


    m_forwardPlan = fftw_plan_dft_3d(m_nPointsTemporal, m_nPointsSpatial, m_nPointsSpatial,
                                     reinterpret_cast<fftw_complex*> (m_real.memptr()),
                                     reinterpret_cast<fftw_complex*> (m_complex.memptr()),
                                     FFTW_FORWARD, FFTW_ESTIMATE);

    m_backwardPlan = fftw_plan_dft_3d(m_nPointsTemporal, m_nPointsSpatial, m_nPointsSpatial,
                                      reinterpret_cast<fftw_complex*> (m_complex.memptr()),
                                      reinterpret_cast<fftw_complex*> (m_real.memptr()),
                                      FFTW_BACKWARD, FFTW_ESTIMATE);

}

Integrator::~Integrator()
{
}


cube Integrator::backwardFFT(cx_cube in)
{
    fftw_execute_dft(m_backwardPlan, reinterpret_cast<fftw_complex*>(in.memptr()),
                     reinterpret_cast<fftw_complex*>(m_real.memptr()));


    //fftShift:
    //shift output to be symmetric around center
    //only in space since temporal should be centered around (0,0)
    for(int i = 0; i < int(m_real.n_slices); i++){
        m_real.slice(i) = FFTHelper::ifftShift(m_real.slice(i));
    }

    m_real *= m_temporalFreqResolution
            * m_spatialFreqResolution
            * m_spatialFreqResolution
            /8./core::pi/core::pi/core::pi;

    return real(m_real);
}


cx_cube Integrator::forwardFFT(cube in)
{
    //fftShift:
    //shift input to be symmetric around (0,0)
    //only in space since temporal is already centered around (0,0)
    for(int i = 0; i < int(m_real.n_slices); i++){
        in.slice(i) = FFTHelper::fftShift(in.slice(i));
    }


    cx_cube tmp = zeros<cx_cube>(in.n_rows, in.n_rows, in.n_slices);
    tmp.set_real(in);

    fftw_execute_dft(m_forwardPlan, reinterpret_cast<fftw_complex*>(tmp.memptr()),
                     reinterpret_cast<fftw_complex*> (m_complex.memptr()));

    return m_complex * m_temporalResolution * m_spatialResolution * m_spatialResolution;
}


//cx_mat Integrator::forwardFFT(mat data) const
//{
//    //fftShift
//    data = FFTHelper::fftShift(data);


//    cx_mat ifftData = zeros<cx_mat>(int(data.n_cols), int(data.n_rows));
//    cx_mat inData = 0*ifftData;
//    inData.set_real(data);
//    int size[2] = {int(data.n_cols), int(data.n_rows)};

//    fftw_complex* in = reinterpret_cast<fftw_complex*> (inData.memptr());
//    fftw_complex* out = reinterpret_cast<fftw_complex*> (ifftData.memptr());
//    fftw_plan plan = fftw_plan_dft(2, size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

//    fftw_execute(plan);
//    fftw_destroy_plan(plan);

//    return ifftData;
//}

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
