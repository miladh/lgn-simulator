#include "gangliondog.h"

GanglionDOG::GanglionDOG(const Config *cfg)
    : Ganglion(cfg)
{

    const Setting & root = cfg->getRoot();
    double A = root["dogSettings"]["A"];
    double a = root["dogSettings"]["a"];
    double B = root["dogSettings"]["B"];
    double b = root["dogSettings"]["b"];

    m_dog = new DOG(A, a, B, b);
    m_tau_rg = root["temporalSettings"]["tau_rg"];
}

GanglionDOG::~GanglionDOG()
{

}

void GanglionDOG::setSpatialImpulseResponse(vec2 r)
{
    m_spatialImpulseResponse = m_dog->real(r);
}

void GanglionDOG::setSpatialImpulseResponseComplex(vec2 k)
{
    m_spatialImpulseResponseComplex = m_dog->complex(k);
}

void GanglionDOG::setTemporalImpulseResponse(double t)
{
    m_temporalImpulseResponse = 1./ m_tau_rg * exp(-t/m_tau_rg)
            * m_math.heaviside(t);
}

void GanglionDOG::setTemporalImpulseResponseComplex(double w)
{
    m_temporalImpulseResponseComplex = 1./ (1 + w*w * m_tau_rg*m_tau_rg);
}

