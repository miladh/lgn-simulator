#include "naturalscene.h"

NaturalScene::NaturalScene(Integrator *integrator, mat scene)
    : Stimulus(integrator)
{
    m_scene = zeros<cx_mat>(m_fourierTransform.slice(0).n_cols,
                            m_fourierTransform.slice(0).n_rows);

    m_sceneFourierTransform = zeros<cx_mat>(m_fourierTransform.slice(0).n_cols,
                                            m_fourierTransform.slice(0).n_rows);
    m_scene.set_real(scene);
}

NaturalScene::~NaturalScene()
{

}

void NaturalScene::computeSpatiotemporal()
{

    for(int k = 0; k < int(m_spatioTemporal.n_slices); k++){
        m_spatioTemporal.slice(k) = real(m_scene);
    }

}

void NaturalScene::computeFourierTransform()
{
    int size[2] = {int(m_scene.n_cols), int(m_scene.n_rows)};

    fftw_complex* in = reinterpret_cast<fftw_complex*> (m_scene.memptr());
    fftw_complex* out = reinterpret_cast<fftw_complex*> (m_sceneFourierTransform.memptr());
    fftw_plan plan = fftw_plan_dft(2, size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_execute(plan);
    fftw_destroy_plan(plan);
//    m_sceneFourierTransform *= m_dk * m_dk;

    for(int k = 0; k < int(m_spatioTemporal.n_slices); k++){
        m_fourierTransform.slice(k) = m_sceneFourierTransform
                * Functions::delta(0, m_temporalFreqs[k])
                /m_integrator->temporalFreqResolution();
    }
}



double NaturalScene::valueAtPoint(vec2 rVec, double t)
{
}

double NaturalScene::fourierTransformAtFrequency(vec2 k, double w)
{
}


NaturalScene createNaturalSceneStimulus(Integrator *integrator, const Config *cfg)
{
    const Setting & root = cfg->getRoot();
    string scenePath = root["stimuliSettings"]["scenePath"];


    cv::Mat cvMat = cv::imread(scenePath, 0);
    if (cvMat.empty())
    {
        cout << "Cannot open image!" << endl;
    }

    cvMat.convertTo(cvMat, CV_64F);
    mat scene(reinterpret_cast<double*>(cvMat.data), cvMat.rows, cvMat.cols);

    return NaturalScene(integrator, fliplr(scene));
}
