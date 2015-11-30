#include "outputmanager.h"
#include <unistd.h>

OutputManager::OutputManager(const Config *cfg)
    : m_cfg(cfg)
{
    const Setting & root = m_cfg->getRoot();
    string outputFilePath = root["fileManagerSettings"]["outputFilePath"];
    m_outputFileName << outputFilePath << "/output.h5";
    m_output = new H5File (m_outputFileName.str(), H5F_ACC_TRUNC);

    initialize();

}
OutputManager::~OutputManager()
{
    delete m_output;
}

void OutputManager::initialize()
{


    Group rootGroup = m_output->openGroup("/");
    const Setting & root = m_cfg->getRoot();


    double dt = root["integratorSettings"]["dt"];
    double ds = root["integratorSettings"]["ds"];
    int nSteps = root["integratorSettings"]["nt"];
    int nPoints = root["integratorSettings"]["ns"];
    nSteps = pow(2, nSteps);
    nPoints = pow(2, nPoints);

    Attribute dt_a(rootGroup.createAttribute("dt",PredType::NATIVE_DOUBLE, H5S_SCALAR));
    Attribute ds_a(rootGroup.createAttribute("ds",PredType::NATIVE_DOUBLE, H5S_SCALAR));
    Attribute nSteps_a(rootGroup.createAttribute("nSteps",PredType::NATIVE_INT, H5S_SCALAR));
    Attribute nPoints_a(rootGroup.createAttribute("nPoints",PredType::NATIVE_INT, H5S_SCALAR));


    dt_a.write(PredType::NATIVE_DOUBLE, &dt);
    ds_a.write(PredType::NATIVE_DOUBLE, &ds);
    nSteps_a.write(PredType::NATIVE_INT, &nSteps);
    nPoints_a.write(PredType::NATIVE_INT, &nPoints);

    m_dataset.reserve(dt);
    m_dataset.reserve(ds);
    m_dataset.reserve(nSteps);
    m_dataset.reserve(nPoints);
}


void OutputManager::writeResponse(const vector<Neuron*> &neurons,
                                  const Stimulus &stimuli){



    // Write stimuli
    Group stim = m_output->createGroup("/stimulus");
    cube realStim = stimuli.spatioTemporal();
    cube complexStim = real(stimuli.fourierTransform());

    writeDataSet(realStim, &stim, "spatioTemporal");
    writeDataSet(complexStim, &stim, "fourierTransform");



    // Write neurons
    for(const Neuron *neuron : neurons){

        cube realResponse = neuron->response();
        cube complexResponse = real(neuron->responseFT());

        cube realImpulseResponse = neuron->impulseResponse();
        cube complexImpulseResponse = real(neuron->impulseResponseFourierTransform());

        string cellGroupName = neuron->cellType();
        Group cellGroup = m_output->createGroup(cellGroupName);


        //write response:
        Group res = m_output->createGroup(cellGroupName+"/response");

        writeDataSet(realResponse, &res, "spatioTemporal");
        writeDataSet(complexResponse, &res, "fourierTransform");

        //write impulse response:
        Group impRes = m_output->createGroup(cellGroupName+"/impulseResponse");
        writeDataSet(realImpulseResponse, &impRes, "spatioTemporal");
        writeDataSet(complexImpulseResponse, &impRes, "fourierTransform");

    }


}

void OutputManager::writeDataSet(const cube data, Group* group, string name)
{
    hsize_t dim[3] = {data.n_slices, data.n_cols, data.n_rows};
    DataSpace space(3, dim);
    DataSet dataset(group->createDataSet(name, PredType::NATIVE_DOUBLE, space));
    dataset.write(data.memptr(), PredType::NATIVE_DOUBLE);

}

