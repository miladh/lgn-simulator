#include "outputmanager.h"
#include <unistd.h>

using namespace edog;


OutputManager::OutputManager(const YAML::Node *cfg)
    : m_cfg(cfg)
{

    string outputFile = (*m_cfg)["outputFile"].as<std::string>();
    m_output = new H5File (outputFile, H5F_ACC_TRUNC);

    initialize();

}
OutputManager::~OutputManager()
{
    delete m_output;
}

void OutputManager::initialize()
{
    Group rootGroup = m_output->openGroup("/");
    double dt = (*m_cfg)["dt"].as<double>();
    int nSteps = (*m_cfg)["nt"].as<int>();
    int nPoints =(*m_cfg)["ns"].as<int>();


    double ds = 1.0/nPoints;
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



void OutputManager::writeResponse(const Neuron* neuron)
{

    cube realResponse = neuron->response();
    cube complexResponse = real(neuron->responseFT());

    string cellGroupName = neuron->cellType();
    herr_t status = H5Eset_auto1(NULL, NULL);
    status = H5Gget_objinfo (m_output->getId(), cellGroupName.c_str(), 0, NULL);

    if (!status == 0){
//        cout << "writeResponse: Cell group doesn't exist.....creating group" << endl;
        Group cellGroup = m_output->createGroup(cellGroupName);
    }

    //write response:
    Group res = m_output->createGroup(cellGroupName+"/response");

    writeDataSet(realResponse, &res, "spatioTemporal");
    writeDataSet(complexResponse, &res, "fourierTransform");

}

void OutputManager::writeImpulseResponse(const Neuron* neuron)
{

    cube realImpulseResponse = neuron->impulseResponse();
    cube complexImpulseResponse = real(neuron->impulseResponseFourierTransform());

    string cellGroupName = neuron->cellType();
    herr_t status = H5Eset_auto1(NULL, NULL);
    status = H5Gget_objinfo (m_output->getId(), cellGroupName.c_str(), 0, NULL);

    if (!status == 0){
//        cout << "writeResponse: Cell group doesn't exist.....creating group" << endl;
        Group cellGroup = m_output->createGroup(cellGroupName);
    }

    //write impulse response:
    Group impRes = m_output->createGroup(cellGroupName+"/impulseResponse");
    writeDataSet(realImpulseResponse, &impRes, "spatioTemporal");
    writeDataSet(complexImpulseResponse, &impRes, "fourierTransform");

}


void OutputManager::writeStimulus(const Stimulus* stimuli)
{

    // Write stimuli
    Group stim = m_output->createGroup("/stimulus");

    if(stimuli->type() == "grating"){
        double C = (*m_cfg)["C"].as<double>();
        string mask = (*m_cfg)["mask"].as<string>();
        double maskSize =(*m_cfg)["maskSize"].as<double>();


        Attribute C_a(stim.createAttribute("C",PredType::NATIVE_DOUBLE, H5S_SCALAR));
        Attribute mask_a(stim.createAttribute("mask", StrType(PredType::C_S1, 64), H5S_SCALAR));
        Attribute maskSize_a(stim.createAttribute("maskSize",PredType::NATIVE_DOUBLE, H5S_SCALAR));


        C_a.write(PredType::NATIVE_DOUBLE, &C);
        mask_a.write( StrType(PredType::C_S1, 64), (&mask)->c_str());
        maskSize_a.write(PredType::NATIVE_DOUBLE, &maskSize);
    }


    cube realStim = stimuli->spatioTemporal();
    cube complexStim = real(stimuli->fourierTransform());
    writeDataSet(realStim, &stim, "spatioTemporal");
    writeDataSet(complexStim, &stim, "fourierTransform");


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


    hsize_t dim[3] = {data.n_slices, data.n_rows, data.n_cols};
    DataSpace space(3, dim);
    DataSet dataset(group->createDataSet(name, PredType::NATIVE_DOUBLE, space));

    dataset.write(data.memptr(), PredType::NATIVE_DOUBLE);

}

