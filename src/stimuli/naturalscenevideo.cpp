#include "naturalscenevideo.h"

using namespace lgnSimulator;


NaturalSceneVideo::NaturalSceneVideo(Integrator* const integrator, string sceneFilename)
    : Stimulus(integrator)
    , m_sceneFilename(sceneFilename)

{
    m_scene = zeros<cube>(m_integrator->nPointsSpatial(),
                            m_integrator->nPointsSpatial(),
                            m_integrator->nPointsTemporal());

    m_type = "naturalSceneVideo";
    readScene();

}

NaturalSceneVideo::~NaturalSceneVideo()
{

}

void NaturalSceneVideo::computeSpatiotemporal()
{
    m_spatiotemporal = real(m_scene);
}

void NaturalSceneVideo::computeFourierTransform()
{
//    cx_mat tmp = m_integrator->forwardFFT(m_scene);
//     m_fourierTransform = m_integrator->forwardFFT(m_scene);
}

void NaturalSceneVideo::readScene()
{

    cv::VideoCapture capture(m_sceneFilename);
    cv::Mat frame;

    if( !capture.isOpened() ){
        throw "Error when reading video";
    }

//    cout  <<  capture.get(CV_CAP_PROP_FRAME_COUNT) << endl;
//    cout <<  capture.get(CV_CAP_PROP_FPS) << endl;

    int i = 0;
    for( ; ; )
    {
        capture.read(frame);
        if(frame.empty() || i > m_integrator->nPointsTemporal() - 1 ){
            break;
        }
        cvtColor(frame, frame, CV_BGR2GRAY);

        //Scaling
        if(frame.rows != m_integrator->nPointsSpatial()
                || frame.cols != m_integrator->nPointsSpatial()){
            cv::Size size(m_integrator->nPointsSpatial(), m_integrator->nPointsSpatial());
            cv::resize(frame, frame, size);
        }

        //Convert to arma mat
        frame.convertTo(frame, CV_64F);
        mat scene(reinterpret_cast<double*>(frame.data), frame.rows, frame.cols);
        m_scene.slice(i).set_real(scene);
        m_scene.slice(i) = fliplr(m_scene.slice(i));
        i+=1;
    }


}



unique_ptr<NaturalSceneVideo> createNaturalSceneVideoStimulus(Integrator* const integrator, const YAML::Node &cfg)
{
    //Read file
    string sceneFilename = cfg["videoFilename"].as<string>();

    return unique_ptr<NaturalSceneVideo>(new NaturalSceneVideo(integrator, sceneFilename));
}
