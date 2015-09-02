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


void OutputManager::writeResponse(const int state,
                                  const vector<Neuron*> &neurons,
                                  const Stimuli &stimuli)
{
    stringstream stateIndex;
    stateIndex << "state" << setw(4) << setfill('0')  << state;
    string stateName = "/"+stateIndex.str();
    Group* group = new Group( m_output->createGroup(stateName));


    // Write stimuli
    Group* stim = new Group( m_output->createGroup(stateName+"/stimuli"));
    mat realStim = stimuli.spatial();
    mat complexStim = real(stimuli.frequency());

    writeDataSet(realStim, stim, "real");
    writeDataSet(complexStim, stim, "complex");


    // Write neurons
    for(const Neuron *neuron : neurons){

        mat realResponse = neuron->response();
        mat complexResponse = neuron->responseComplex();

        mat realImpulseResponse = neuron->impulseResponse();
        mat complexImpulseResponse = neuron->impulseResponseComplex();

        string cellGroupName = stateName+"/"+neuron->cellType();
        Group* cellGroup = new Group( m_output->createGroup(cellGroupName));

        //write response:
        Group* res = new Group( m_output->createGroup(cellGroupName+"/response"));

        writeDataSet(realResponse, res, "real");
        writeDataSet(complexResponse, res, "complex");

        //write impulse response:
        Group* impRes =
                new Group( m_output->createGroup(cellGroupName+"/impulseResponse"));
        writeDataSet(realImpulseResponse, impRes, "real");
        writeDataSet(complexImpulseResponse, impRes, "complex");



    }

}


void OutputManager::writeDataSet(const mat data, Group* group, string name)
{
    hsize_t dim[2] = {data.n_cols, data.n_rows};
    DataSpace space(2, dim);
    DataSet dataset(group->createDataSet(name, PredType::NATIVE_DOUBLE, space));
    dataset.write(data.memptr(), PredType::NATIVE_DOUBLE);

}

