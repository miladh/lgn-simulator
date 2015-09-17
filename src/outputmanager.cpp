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

    int nSteps = root["dynamicSettings"]["nSteps"];
    double dt = root["dynamicSettings"]["dt"];

    Attribute nSteps_a(rootGroup.createAttribute("nSteps", PredType::NATIVE_INT, H5S_SCALAR));
    nSteps_a.write(PredType::NATIVE_INT, &nSteps);

    Attribute dt_a(rootGroup.createAttribute("dt", PredType::NATIVE_DOUBLE, H5S_SCALAR));
    dt_a.write(PredType::NATIVE_DOUBLE, &dt);

    m_dataset.reserve(nSteps);
}


void OutputManager::writeResponse(const vector<Neuron*> &neurons,
                                  const Stimuli &stimuli){



    // Write stimuli
    Group stim = m_output->createGroup("/stimuli");
    cube realStim = stimuli.spatioTemporal();
    cube complexStim = real(stimuli.fourierTransform());

    //   cout <<realStim.size() << endl;
    writeDataSet(realStim, &stim, "real");
    writeDataSet(complexStim, &stim, "complex");



    // Write neurons
    for(const Neuron *neuron : neurons){

        cube realResponse = neuron->response();
        cube complexResponse = real(neuron->responseFT());

        cube realImpulseResponse = neuron->impulseResponse();
        cube complexImpulseResponse = real(neuron->impulseResponseFT());

        string cellGroupName = neuron->cellType();
        Group cellGroup = m_output->createGroup(cellGroupName);

        //write response:
        Group res = m_output->createGroup(cellGroupName+"/response");

        writeDataSet(realResponse, &res, "real");
        writeDataSet(complexResponse, &res, "complex");

        //write impulse response:
        Group impRes = m_output->createGroup(cellGroupName+"/impulseResponse");
        writeDataSet(realImpulseResponse, &impRes, "real");
        writeDataSet(complexImpulseResponse, &impRes, "complex");

    }


}

void OutputManager::writeDataSet(const cube data, Group* group, string name)
{
    hsize_t dim[3] = {data.n_slices, data.n_cols, data.n_rows};
    DataSpace space(3, dim);
    DataSet dataset(group->createDataSet(name, PredType::NATIVE_DOUBLE, space));
    dataset.write(data.memptr(), PredType::NATIVE_DOUBLE);

}

