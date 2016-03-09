#include "outputmanager.h"
#include <unistd.h>

#include "stimuli/grating/grating.h"
#include "stimuli/naturalscene.h"
#include "stimuli/naturalscenevideo.h"

using namespace lgnSimulator;


OutputManager::OutputManager(const string& filename)
{
    m_output = new H5File (filename, H5F_ACC_TRUNC);
}
OutputManager::~OutputManager()
{
    delete m_output;
}

void OutputManager::writeIntegratorProperties(const Integrator &integrator)
{
    // Write stimuli
    Group integratorGroup = m_output->createGroup("/integrator");
    vector <DataSet *>dataset;

    int Nt = integrator.nPointsTemporal();
    int Ns = integrator.nPointsSpatial();
    double dt = integrator.temporalResolution();
    double ds = integrator.spatialResolution();
    double dw = integrator.temporalFreqResolution();
    double dk = integrator.spatialFreqResolution();
    double ws = integrator.temporalSamplingFreq();
    double ks = integrator.spatialSamplingFreq();
    double T = integrator.timeInterval();
    double L = integrator.lengthInterval();

    vec timeVec = integrator.timeVec();
    vec spatialVec = integrator.spatialVec();
    vec temporalFreqVec = integrator.temporalFreqVec();
    vec spatialFreqVec = integrator.spatialFreqVec();


    Attribute Nt_a(integratorGroup.createAttribute("nPointsTemporal",PredType::NATIVE_INT, H5S_SCALAR));
    Attribute Ns_a(integratorGroup.createAttribute("nPointsSpatial",PredType::NATIVE_INT, H5S_SCALAR));
    Attribute dt_a(integratorGroup.createAttribute("temporalResolution",PredType::NATIVE_DOUBLE, H5S_SCALAR));
    Attribute ds_a(integratorGroup.createAttribute("spatialResolution",PredType::NATIVE_DOUBLE, H5S_SCALAR));
    Attribute dw_a(integratorGroup.createAttribute("temporalFreqResolution",PredType::NATIVE_DOUBLE, H5S_SCALAR));
    Attribute dk_a(integratorGroup.createAttribute("spatialFreqResolution",PredType::NATIVE_DOUBLE, H5S_SCALAR));
    Attribute ws_a(integratorGroup.createAttribute("temporalSamplingFreq",PredType::NATIVE_DOUBLE, H5S_SCALAR));
    Attribute ks_a(integratorGroup.createAttribute("spatialSamplingFreq",PredType::NATIVE_DOUBLE, H5S_SCALAR));
    Attribute T_a(integratorGroup.createAttribute("timeInterval",PredType::NATIVE_DOUBLE, H5S_SCALAR));
    Attribute L_a(integratorGroup.createAttribute("lengthInterval",PredType::NATIVE_DOUBLE, H5S_SCALAR));

    Nt_a.write(PredType::NATIVE_INT, &Nt);
    Ns_a.write(PredType::NATIVE_INT, &Ns);
    dt_a.write(PredType::NATIVE_DOUBLE, &dt);
    ds_a.write(PredType::NATIVE_DOUBLE, &ds);
    dw_a.write(PredType::NATIVE_DOUBLE, &dw);
    dk_a.write(PredType::NATIVE_DOUBLE, &dk);
    ws_a.write(PredType::NATIVE_DOUBLE, &ws);
    ks_a.write(PredType::NATIVE_DOUBLE, &ks);
    T_a.write(PredType::NATIVE_DOUBLE, &T);
    L_a.write(PredType::NATIVE_DOUBLE, &L);

    dataset.reserve(Nt);
    dataset.reserve(Ns);
    dataset.reserve(dt);
    dataset.reserve(ds);
    dataset.reserve(dw);
    dataset.reserve(dk);
    dataset.reserve(ws);
    dataset.reserve(ks);
    dataset.reserve(T);
    dataset.reserve(L);

    hsize_t dim_t[1] = {timeVec.n_elem};
    hsize_t dim_s[1] = {spatialVec.n_elem};

    DataSpace space_t(1, dim_t);
    DataSpace space_s(1, dim_s);

    DataSet timeVec_ds(integratorGroup.createDataSet("timeVec", PredType::NATIVE_DOUBLE, space_t));
    DataSet temporalFreqVec_ds(integratorGroup.createDataSet("temporalFreqVec", PredType::NATIVE_DOUBLE, space_t));

    DataSet spatialVec_ds(integratorGroup.createDataSet("spatialVec", PredType::NATIVE_DOUBLE, space_s));
    DataSet spatialFreqVec_ds(integratorGroup.createDataSet("spatialFreqVec", PredType::NATIVE_DOUBLE, space_s));

    timeVec_ds.write(timeVec.memptr(), PredType::NATIVE_DOUBLE);
    temporalFreqVec_ds.write(temporalFreqVec.memptr(), PredType::NATIVE_DOUBLE);
    spatialVec_ds.write(spatialVec.memptr(), PredType::NATIVE_DOUBLE);
    spatialFreqVec_ds.write(spatialFreqVec.memptr(), PredType::NATIVE_DOUBLE);
}



void OutputManager::writeResponse(const Neuron& neuron)
{

    cube response = neuron.response();
    cube responseFT = real(neuron.responseFT());

    string cellGroupName = neuron.type();
    herr_t status = H5Eset_auto1(NULL, NULL);
    status = H5Gget_objinfo (m_output->getId(), cellGroupName.c_str(), 0, NULL);

    if (!status == 0){
        Group cellGroup = m_output->createGroup(cellGroupName);
    }

    //write response:
    Group res = m_output->createGroup(cellGroupName+"/response");

    writeDataSet(response, &res, "spatioTemporal");
    writeDataSet(responseFT, &res, "fourierTransform");

}

void OutputManager::writeImpulseResponse(const Neuron& neuron)
{

    cube impulseResponse = neuron.impulseResponse();
    cube impulseResponseFT = real(neuron.impulseResponseFourierTransform());

    string cellGroupName = neuron.type();
    herr_t status = H5Eset_auto1(NULL, NULL);
    status = H5Gget_objinfo (m_output->getId(), cellGroupName.c_str(), 0, NULL);

    if (!status == 0){
        Group cellGroup = m_output->createGroup(cellGroupName);
    }

    //write impulse response:
    Group impRes = m_output->createGroup(cellGroupName+"/impulseResponse");
    writeDataSet(impulseResponse, &impRes, "spatioTemporal");
    writeDataSet(impulseResponseFT, &impRes, "fourierTransform");

}


void OutputManager::writeStimulus(const Stimulus* stimulus)
{

    // Write stimuli
    Group stim = m_output->createGroup("/stimulus");
    string type = stimulus->type();

    Attribute type_a(stim.createAttribute("type", StrType(PredType::C_S1, 64), H5S_SCALAR));
    type_a.write( StrType(PredType::C_S1, 64), (&type)->c_str());


    if (const  Grating * gratingStimulus = dynamic_cast<const Grating*>(stimulus) ) {
        double C = gratingStimulus->contrast();
        double maskSize = gratingStimulus->maskSize();
        string mask = gratingStimulus->mask();
        vec2 k = gratingStimulus->k();
        double w = gratingStimulus->w();


        Attribute C_a(stim.createAttribute("C",PredType::NATIVE_DOUBLE, H5S_SCALAR));
        Attribute mask_a(stim.createAttribute("mask", StrType(PredType::C_S1, 64), H5S_SCALAR));
        Attribute maskSize_a(stim.createAttribute("maskSize",PredType::NATIVE_DOUBLE, H5S_SCALAR));
        Attribute kx_a(stim.createAttribute("kx",PredType::NATIVE_DOUBLE, H5S_SCALAR));
        Attribute ky_a(stim.createAttribute("ky",PredType::NATIVE_DOUBLE, H5S_SCALAR));
        Attribute w_a(stim.createAttribute("w",PredType::NATIVE_DOUBLE, H5S_SCALAR));



        C_a.write(PredType::NATIVE_DOUBLE, &C);
        mask_a.write( StrType(PredType::C_S1, 64), (&mask)->c_str());
        maskSize_a.write(PredType::NATIVE_DOUBLE, &maskSize);
        kx_a.write(PredType::NATIVE_DOUBLE, &k(0));
        ky_a.write(PredType::NATIVE_DOUBLE, &k(1));
        w_a.write(PredType::NATIVE_DOUBLE, &w);
    }


    cube realStim = stimulus->spatioTemporal();
    cube complexStim = real(stimulus->fourierTransform());
    writeDataSet(realStim, &stim, "spatioTemporal");
    writeDataSet(complexStim, &stim, "fourierTransform");

}





void OutputManager::writeDataSet(const cube data, Group* group, string name)
{

    hsize_t dim[3] = {data.n_slices, data.n_rows, data.n_cols};
    DataSpace space(3, dim);
    DataSet dataset(group->createDataSet(name, PredType::NATIVE_DOUBLE, space));

    dataset.write(data.memptr(), PredType::NATIVE_DOUBLE);

}

