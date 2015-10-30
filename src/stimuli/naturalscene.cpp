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
        m_spatioTemporal.slice(k) = real(m_scene)/**cos(m_integrator->temporalFreqVec()[1])*/;
    }

}

void NaturalScene::computeFourierTransform()
{
    m_sceneFourierTransform = m_integrator->forwardFFT(m_scene);

    for(int k = 0; k < int(m_spatioTemporal.n_slices); k++){
        m_fourierTransform.slice(k) = m_sceneFourierTransform
                * Functions::delta(0, m_temporalFreqs[k])
                /m_integrator->temporalFreqResolution();
    }
}


NaturalScene createNaturalSceneStimulus(Integrator *integrator, const Config *cfg)
{
    //Read file
    const Setting & root = cfg->getRoot();
    string scenePath = root["stimuliSettings"]["scenePath"];
    int ns = root["integratorSettings"]["ns"];
    int nt = root["integratorSettings"]["nt"];

    int Ns = pow(2,ns);
    int Nt = pow(2,nt);


    cv::Mat cvMat = cv::imread(scenePath, 0);
    if (cvMat.empty())
    {
        cout << "Cannot open image!" << endl;
    }

    //    cv::imshow("image", cvMat);
    //    cv::waitKey(0);


    //Scaling
    if(cvMat.rows != Ns ||cvMat.cols != Ns){
        cv::Size size(Ns, Ns);
        cv::resize(cvMat, cvMat, size);
    }


    //Convert to arma mat
    cvMat.convertTo(cvMat, CV_64F);
    mat scene(reinterpret_cast<double*>(cvMat.data), cvMat.rows, cvMat.cols);

    return NaturalScene(integrator, scene);
}
