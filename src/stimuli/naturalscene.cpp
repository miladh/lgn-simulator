#include "naturalscene.h"

using namespace edog;


NaturalScene::NaturalScene(Integrator *integrator, string sceneFilename)
    : Stimulus(integrator)
    , m_sceneFilename(sceneFilename)

{
    m_scene = zeros<cx_mat>(m_integrator->nPointsSpatial(), m_integrator->nPointsSpatial());
    m_sceneFourierTransform = zeros<cx_mat>(m_integrator->nPointsSpatial(),
                                            m_integrator->nPointsSpatial());
    m_type =  "naturalScene";
    readScene();
}

NaturalScene::~NaturalScene()
{

}

void NaturalScene::computeSpatiotemporal()
{

    for(int k = 0; k < int(m_spatioTemporal.n_slices); k++){
        m_spatioTemporal.slice(k) = real(m_scene)
                * temporalValueAtPoint(m_timeVec[k]);
    }

}

void NaturalScene::computeFourierTransform()
{
    m_sceneFourierTransform = m_integrator->forwardFFT(m_scene);

    for(int k = 0; k < int(m_spatioTemporal.n_slices); k++){
        m_fourierTransform.slice(k) = m_sceneFourierTransform
                * fourierTransformAtTemporalFrequency(m_temporalFreqs[k]);
    }
}

void NaturalScene::readScene()
{
    cv::Mat cvMat = cv::imread(m_sceneFilename, 0);
    if (cvMat.empty())
    {
        throw "Cannot open image!";
    }

    //    cv::imshow("image", cvMat);
    //    cv::waitKey(0);

    //Scaling
    if(cvMat.rows != m_integrator->nPointsSpatial()
            || cvMat.cols != m_integrator->nPointsSpatial()){
        cv::Size size(m_integrator->nPointsSpatial(), m_integrator->nPointsSpatial());
        cv::resize(cvMat, cvMat, size);
    }


    //Convert to arma mat
    cvMat.convertTo(cvMat, CV_64F);
    mat scene(reinterpret_cast<double*>(cvMat.data), cvMat.rows, cvMat.cols);
    m_scene.set_real(scene);
}
