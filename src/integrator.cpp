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


    //fft plans
    m_out = zeros<cx_cube>(m_nPointsSpatial, m_nPointsSpatial, m_nPointsTemporal);
    m_in  = zeros<cx_cube>(m_nPointsSpatial, m_nPointsSpatial, m_nPointsTemporal);
    m_in_2d = zeros<cx_mat>(m_nPointsSpatial, m_nPointsSpatial);
    m_out_2d = zeros<cx_mat>(m_nPointsSpatial, m_nPointsSpatial);


    m_backwardPlan = fftw_plan_dft_3d(m_nPointsTemporal, m_nPointsSpatial, m_nPointsSpatial,
                                      reinterpret_cast<fftw_complex*> (m_out.memptr()),
                                      reinterpret_cast<fftw_complex*> (m_in.memptr()),
                                      FFTW_BACKWARD, FFTW_MEASURE);

    m_forwardPlan = fftw_plan_dft_3d(m_nPointsTemporal, m_nPointsSpatial, m_nPointsSpatial,
                                     reinterpret_cast<fftw_complex*> (m_in.memptr()),
                                     reinterpret_cast<fftw_complex*> (m_out.memptr()),
                                     FFTW_FORWARD, FFTW_MEASURE);


    m_forwardPlan_2d = fftw_plan_dft_2d(m_nPointsSpatial, m_nPointsSpatial,
                                     reinterpret_cast<fftw_complex*> (m_in_2d.memptr()),
                                     reinterpret_cast<fftw_complex*> (m_out_2d.memptr()),
                                     FFTW_FORWARD, FFTW_MEASURE);
}

Integrator::~Integrator()
{
}


cube Integrator::backwardFFT(const cx_cube& in)
{
    m_in = in;
    fftw_execute_dft(m_backwardPlan, reinterpret_cast<fftw_complex*>(m_in.memptr()),
                     reinterpret_cast<fftw_complex*>(m_out.memptr()));

    //fftShift:
    //shift output to be symmetric around center
    //only in space since temporal should be centered around (0,0)
    for(int i = 0; i < int(m_out.n_slices); i++){
        m_out.slice(i) = FFTHelper::ifftShift(m_out.slice(i));
    }

    m_out *= m_temporalFreqResolution
            * m_spatialFreqResolution
            * m_spatialFreqResolution
            /8./core::pi/core::pi/core::pi;

    return real(m_out);
}


cx_cube Integrator::forwardFFT(const cube &in)
{

    //fftShift:
    //shift input to be symmetric around (0,0)
    //only in space since temporal is already centered around (0,0)
    for(int i = 0; i < int(m_in.n_slices); i++){
        m_in.slice(i).set_real(FFTHelper::fftShift(in.slice(i)));
        m_in.slice(i).set_imag(zeros(m_in.n_rows, m_in.n_cols));
    }

    fftw_execute_dft(m_forwardPlan, reinterpret_cast<fftw_complex*>(m_in.memptr()),
                     reinterpret_cast<fftw_complex*> (m_out.memptr()));

    return m_out * m_temporalResolution * m_spatialResolution * m_spatialResolution;
}


cx_mat Integrator::forwardFFT(const mat &in)
{
    //fftShift:
    //shift input to be symmetric around (0,0)
    m_in_2d.set_real(FFTHelper::fftShift(in));
    m_in_2d.set_imag(zeros(in.n_rows, in.n_cols));

    fftw_execute_dft(m_forwardPlan_2d, reinterpret_cast<fftw_complex*>(m_in_2d.memptr()),
                     reinterpret_cast<fftw_complex*> (m_out_2d.memptr()));

    return m_out_2d * m_spatialResolution * m_spatialResolution;;
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
