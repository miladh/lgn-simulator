#include "fullfieldgrating.h"

using namespace lgnSimulator;


FullFieldGrating::FullFieldGrating(const Integrator &integrator,
                                   vec2 kd, double wd, double contrast)
    : Grating(integrator, kd, wd, contrast, 0)
{
    m_mask = "none";
    m_peak = 1.0
            / m_integrator.temporalFreqResolution()
            / m_integrator.spatialFreqResolution()
            / m_integrator.spatialFreqResolution();

}

FullFieldGrating::~FullFieldGrating()
{
}

double FullFieldGrating::valueAtPoint(vec2 rVec, double t) const
{
//    double phi = 2*core::pi * t/m_integrator.timeInterval();
//    vec2 k_rot = 0*m_k;
//    k_rot(0) = cos(phi)*m_k(0) - sin(phi)*m_k(1);
//    k_rot(1) = cos(phi)*m_k(1) + sin(phi)*m_k(0);


//    double phi = core::pi/3;
//    double k_mag = sqrt(dot(m_k, m_k));
//    vec2 k_rot = {k_mag*cos(phi), k_mag*sin(phi)};
//    double s = m_contrast * cos(dot(k_rot, rVec) - m_w * t);

        double s = m_contrast * cos(dot(m_k, rVec) - m_w * t );
    return s;
}

complex<double> FullFieldGrating::fourierTransformAtFrequency(vec2 k, double w) const
{
//    double phi = core::pi/3;
//    double k_mag = sqrt(dot(m_k, m_k));
//    vec2 k_rot = {k_mag*cos(phi), k_mag*sin(phi)};

//    k_rot(0) = Special::nearestValue(m_spatialFreqs, k_rot(0));
//    k_rot(1) = Special::nearestValue(m_spatialFreqs, k_rot(1));

//    cout << k_rot.t() << endl;

//    double s = (Special::delta(k, k_rot) * Special::delta(w,m_w)
//                + Special::delta(k, -k_rot) * Special::delta(w, -m_w));



    double s = (Special::delta(k, m_k) * Special::delta(w,m_w)
                + Special::delta(k, -m_k) * Special::delta(w, -m_w));

    return 4.*core::pi*core::pi*core::pi * m_contrast * m_peak * s;
}


