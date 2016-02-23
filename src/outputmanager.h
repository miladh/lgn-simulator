#ifndef OUTPUTMANAGER_H
#define OUTPUTMANAGER_H

#include <iostream>
#include <armadillo>
#include <H5Cpp.h>

#include "stimuli/stimulus.h"
#include "neurons/neuron.h"


using namespace arma;
using namespace std;
using namespace H5;

namespace lgnSimulator {
class OutputManager
{
public:
    OutputManager(const string &filename);
    ~OutputManager();

    void writeIntegratorProperties(const Integrator &integrator);
    void writeResponse(const Neuron &neuron);
    void writeImpulseResponse(const Neuron &neuron);
    void writeStimulus(const Stimulus *stimuli);

private:
    H5File *m_output;
    void writeDataSet(const cube dataset, Group* group, string name);



};
}
#endif // OUTPUTMANAGER_H
