#ifndef OUTPUTMANAGER_H
#define OUTPUTMANAGER_H

#include <iostream>
#include <armadillo>
#include <H5Cpp.h>
#include <boost/filesystem.hpp>

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

    void copyConfigFile(const string &configFilename);

    void writeIntegratorProperties(const Integrator &integrator);
    void writeStimulusProperties(const Stimulus *stimulus);

    void writeStimulus(const Stimulus *stimuli);
    void writeStimulusFourierTransform(const Stimulus *stimuli);

    void writeResponse(const Neuron &neuron);
    void writeResponseFourierTransform(const Neuron &neuron);

    void writeImpulseResponse(const Neuron &neuron);
    void writeImpulseResponseFourierTransform(const Neuron &neuron);


private:
    H5File *m_output;
    void writeDataSet(const fcube dataset, Group* group, string name);
    Group createGroupIfNotExist(const string groupName);


};
}
#endif // OUTPUTMANAGER_H
