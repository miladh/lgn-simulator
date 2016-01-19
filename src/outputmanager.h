#ifndef OUTPUTMANAGER_H
#define OUTPUTMANAGER_H

#include <iostream>
#include <armadillo>
#include <H5Cpp.h>
#include <yaml-cpp/yaml.h>

#include "stimuli/stimuli.h"
#include "neurons/neuron.h"


using namespace arma;
using namespace std;
using namespace H5;

namespace edog {
class OutputManager
{
public:
    OutputManager(const YAML::Node *cfg);
    ~OutputManager();

    void writeResponse(const vector<Neuron *> &neurons,
                       const Stimulus &stimuli);

    void writeResponse(const Neuron *neuron);
    void writeImpulseResponse(const Neuron* neuron);
    void writeStimulus(const Stimulus *stimuli);

private:
    const YAML::Node* m_cfg;

    H5File *m_output;
    stringstream m_outputFileName;
    vector <DataSet *>m_dataset;

    void writeDataSet(const cube dataset, Group* group, string name);


    void initialize();

};
}
#endif // OUTPUTMANAGER_H
