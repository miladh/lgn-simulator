#include "outputmanager.h"
#include <unistd.h>
OutputManager::OutputManager(const Config *cfg)
    : m_cfg(cfg)
{
    //    const Setting & root = m_cfg->getRoot();
    //    string output = root["fileManagerSettings"]["outputFilePath"];
    string outputFilePath = ".";
    m_outputFileName << outputFilePath << "/output.h5";
    m_output = new H5File (m_outputFileName.str(), H5F_ACC_TRUNC);

    initialize();

}
OutputManager::~OutputManager()
{

}

void OutputManager::initialize()
{
    //    const Setting & root = m_cfg->getRoot();
    //    int    nSteps = root["dynamicSettings"]["nSteps"];
    int nSteps = 3.;

    Group rootGroup = m_output->openGroup("/");
    //    Attribute nAtoms_a(rootGroup.createAttribute("nAtoms", PredType::NATIVE_INT, H5S_SCALAR));
    //    nAtoms_a.write(PredType::NATIVE_INT, &m_nAtoms);



    m_dataset.reserve(nSteps);


}


void OutputManager::writeResponse(const int state, const Response &response)
{

    mat realResponse = response.real();
    mat complexResponse = response.complex();

    stringstream stateIndex;
    stateIndex << "state" << setw(4) << setfill('0')  << state;
    Group* group = new Group( m_output->createGroup( "/"+stateIndex.str()));



    hsize_t dim[2] = {realResponse.n_cols, realResponse.n_rows};
    DataSpace space(2, dim);
    DataSet dataset(group->createDataSet("real", PredType::NATIVE_DOUBLE, space));
    dataset.write(realResponse.memptr(), PredType::NATIVE_DOUBLE);

    hsize_t dim1[2] = {complexResponse.n_cols, complexResponse.n_rows};
    DataSpace space1(2, dim1);
    DataSet dataset1(group->createDataSet("complex", PredType::NATIVE_DOUBLE, space1));
    dataset1.write(complexResponse.memptr(), PredType::NATIVE_DOUBLE);


}

