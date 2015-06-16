#ifndef OUTPUTMANAGER_H
#define OUTPUTMANAGER_H

#include <iostream>
#include <armadillo>
#include <H5Cpp.h>
#include <libconfig.h++>

#include "response.h"

using namespace arma;
using namespace std;
using namespace H5;
using namespace libconfig;


class OutputManager
{
public:
    OutputManager(const Config *cfg);
    ~OutputManager();

    void writeResponse(const int state, const Response &response);

private:
    const Config* m_cfg;

    H5File *m_output;
    stringstream m_outputFileName;
    vector <DataSet *>m_dataset;


    void initialize();

};

#endif // OUTPUTMANAGER_H
