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

    void writeIntegratorProperties(const Integrator &integrator);

    void writeResponse(const Neuron &neuron, const bool fourierTransform = true);

    void writeImpulseResponse(const Neuron &neuron,const bool fourierTransform = true);

    void writeStimulusProperties(const Stimulus *stimulus);
    void writeStimulus(const Stimulus *stimuli, const bool fourierTransform = true);

    void copyConfigFile(const string &configFilename);

private:
    H5File *m_output;
    void writeDataSet(const fcube dataset, Group* group, string name);



};
}
#endif // OUTPUTMANAGER_H
