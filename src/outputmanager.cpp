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

    vec realGrid = zeros<vec>(3);
    vec complexGrid = zeros<vec>(3);
    vec domain = zeros<vec>(3);

    int nSteps = root["dynamicSettings"]["nSteps"];
    double dt = root["dynamicSettings"]["dt"];
    const Setting &real = root["gridSettings"]["realGrid"];
    const Setting &complex = root["gridSettings"]["complexGrid"];
    const Setting &integrationDomain = root["gridSettings"]["integrationDomain"];

    for(int i =0; i < 3; i++){
        realGrid[i] = real[i];
        complexGrid[i] = complex[i];
        domain[i] = integrationDomain[i];
    }

    Attribute nSteps_a(rootGroup.createAttribute("nSteps", PredType::NATIVE_INT, H5S_SCALAR));
    nSteps_a.write(PredType::NATIVE_INT, &nSteps);

    Attribute dt_a(rootGroup.createAttribute("dt", PredType::NATIVE_DOUBLE, H5S_SCALAR));
    dt_a.write(PredType::NATIVE_DOUBLE, &dt);


    Attribute rMin(rootGroup.createAttribute("rMin", PredType::NATIVE_DOUBLE, H5S_SCALAR));
    Attribute rMax(rootGroup.createAttribute("rMax", PredType::NATIVE_DOUBLE, H5S_SCALAR));
    Attribute rPoints(rootGroup.createAttribute("rPoints", PredType::NATIVE_DOUBLE, H5S_SCALAR));

    rMin.write(PredType::NATIVE_DOUBLE, &realGrid[0]);
    rMax.write(PredType::NATIVE_DOUBLE, &realGrid[1]);
    rPoints.write(PredType::NATIVE_DOUBLE, &realGrid[2]);


    Attribute kMin(rootGroup.createAttribute("kMin", PredType::NATIVE_DOUBLE, H5S_SCALAR));
    Attribute kMax(rootGroup.createAttribute("kMax", PredType::NATIVE_DOUBLE, H5S_SCALAR));
    Attribute kPoints(rootGroup.createAttribute("kPoints", PredType::NATIVE_DOUBLE, H5S_SCALAR));

    kMin.write(PredType::NATIVE_DOUBLE, &complexGrid[0]);
    kMax.write(PredType::NATIVE_DOUBLE, &complexGrid[1]);
    kPoints.write(PredType::NATIVE_DOUBLE, &complexGrid[2]);


    Attribute min(rootGroup.createAttribute("min", PredType::NATIVE_DOUBLE, H5S_SCALAR));
    Attribute max(rootGroup.createAttribute("max", PredType::NATIVE_DOUBLE, H5S_SCALAR));
    Attribute points(rootGroup.createAttribute("points", PredType::NATIVE_DOUBLE, H5S_SCALAR));

    min.write(PredType::NATIVE_DOUBLE, &domain[0]);
    max.write(PredType::NATIVE_DOUBLE, &domain[1]);
    points.write(PredType::NATIVE_DOUBLE, &domain[2]);


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

