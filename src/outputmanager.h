#ifndef OUTPUTMANAGER_H
#define OUTPUTMANAGER_H

#include <iostream>
#include <armadillo>
#include <H5Cpp.h>
#include <libconfig.h++>

#include "stimuli/stimuli.h"
#include "neurons/neuron.h"

using namespace arma;
using namespace std;
using namespace H5;
using namespace libconfig;


class OutputManager
{
public:
    OutputManager(const Config *cfg);
    ~OutputManager();

    void writeResponse(const vector<Neuron *> &neurons,
                       const Stimulus &stimuli);


//    void writeResponse(const vector<Neuron *> &neurons,
//                       const Stimuli &stimuli);

private:
    const Config* m_cfg;

    H5File *m_output;
    stringstream m_outputFileName;
    vector <DataSet *>m_dataset;

    void writeDataSet(const cube dataset, Group* group, string name);


    void initialize();

};

#endif // OUTPUTMANAGER_H
