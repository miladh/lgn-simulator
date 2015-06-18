#ifndef OUTPUTMANAGER_H
#define OUTPUTMANAGER_H

#include <iostream>
#include <armadillo>
#include <H5Cpp.h>
#include <libconfig.h++>

#include "response.h"
#include "stimuli.h"
#include"impulseResponse.h"

using namespace arma;
using namespace std;
using namespace H5;
using namespace libconfig;


class OutputManager
{
public:
    OutputManager(const Config *cfg);
    ~OutputManager();

    void writeResponse(const int state, const Response &response,
                       const ImpulseResponse &impulseResponse,
                       const Stimuli &stimuli);

private:
    const Config* m_cfg;

    H5File *m_output;
    stringstream m_outputFileName;
    vector <DataSet *>m_dataset;

    void writeDataSet(const mat dataset, Group* group, string name);


    void initialize();

};

#endif // OUTPUTMANAGER_H
