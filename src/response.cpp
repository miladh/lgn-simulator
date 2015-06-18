#include "response.h"

Response::Response(ImpulseResponse *impResFunc, Stimuli *stim, Integrator *integrator,
                   vec3 mesh, vec3 integrationDomain)
    : m_impResFunc(impResFunc)
    , m_stim(stim)
    , m_Integrator(integrator)
    , m_mesh(linspace(mesh[0], mesh[1], mesh[2]))
    , m_domain(integrationDomain)
    , m_response(ones(mesh[2], mesh[2]))
    , m_responseComplex(zeros(mesh[2], mesh[2]))

{

}

Response::~Response()
{

}

void Response::computeComplex(double w)
{

    for(int i = 0; i < int(m_mesh.n_elem); i++){
        for(int j = 0; j < int(m_mesh.n_elem); j++){
            double G = m_impResFunc->edogComplex(m_mesh[i], m_mesh[j], w);
            double s = m_stim->patchGratingComplex(m_mesh[i], m_mesh[j], w);

            m_responseComplex(i,j) = G*s;
        }
    }

}

void Response::compute(double t)
{

    mat stim = 0*m_response;
    mat impRes = 0*m_response;

    m_response = 0*m_response;
    double *w = new double [int(m_domain[2])];
    double *x = new double [int(m_domain[2])];
    gauleg(m_domain[0], m_domain[1], x, w, m_domain[2]);

    for(int i = 0; i < int(m_mesh.n_elem); i++){
        for(int j = 0; j < int(m_mesh.n_elem); j++){

            stim(i,j) = m_stim->patchGrating(m_mesh[i], m_mesh[j], t);
            for(int m = 0; m < int(m_domain[2]); m++){
                for(int n = 0; n < int(m_domain[2]); n++){

                    double G = m_impResFunc->edogComplex(x[m], x[n], m_stim->w() );
                    double s = m_stim->patchGratingComplex(x[m], x[n], m_stim->w() );

                    double value = G * w[m] * w[n] *
                            cos(m_mesh[i]*x[m]+ m_mesh[j]*x[n] - m_stim->w() * t);
                    impRes(i, j) +=  value;
                    m_response(i,j) += value * s;
                }
            }
        }
    }

    m_impResFunc->setReal(impRes);
    m_stim->setReal(stim);


}

mat Response::real() const
{
    return m_response;
}

mat Response::complex() const
{
    return m_responseComplex;
}

