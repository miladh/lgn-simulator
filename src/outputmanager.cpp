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

void OutputManager::copyConfigFile(const string& configFilename)
{

    boost::filesystem::path cfgPath(configFilename);
    boost::filesystem::path p(m_output->getFileName());
    boost::filesystem::path outputDir = p.parent_path();

    string copiedFilePath;
    copiedFilePath  = outputDir.string() + "/_" + cfgPath.filename().string();

    boost::filesystem::copy_file(boost::filesystem::path(configFilename),
                                 boost::filesystem::path(copiedFilePath),
                                 boost::filesystem::copy_option::overwrite_if_exists);

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

    fvec timeVec = conv_to<fvec>::from(integrator.timeVec());
    fvec spatialVec = conv_to<fvec>::from(integrator.spatialVec());
    fvec temporalFreqVec = conv_to<fvec>::from(integrator.temporalFreqVec());
    fvec spatialFreqVec = conv_to<fvec>::from(integrator.spatialFreqVec());


    Attribute Nt_a(integratorGroup.createAttribute("Nt",PredType::NATIVE_INT, H5S_SCALAR));
    Attribute Ns_a(integratorGroup.createAttribute("Ns",PredType::NATIVE_INT, H5S_SCALAR));
    Attribute dt_a(integratorGroup.createAttribute("dt",PredType::NATIVE_DOUBLE, H5S_SCALAR));
    Attribute ds_a(integratorGroup.createAttribute("ds",PredType::NATIVE_DOUBLE, H5S_SCALAR));
    Attribute dw_a(integratorGroup.createAttribute("dw",PredType::NATIVE_DOUBLE, H5S_SCALAR));
    Attribute dk_a(integratorGroup.createAttribute("dk",PredType::NATIVE_DOUBLE, H5S_SCALAR));
    Attribute ws_a(integratorGroup.createAttribute("ws",PredType::NATIVE_DOUBLE, H5S_SCALAR));
    Attribute ks_a(integratorGroup.createAttribute("ks",PredType::NATIVE_DOUBLE, H5S_SCALAR));
    Attribute T_a(integratorGroup.createAttribute("T",PredType::NATIVE_DOUBLE, H5S_SCALAR));
    Attribute L_a(integratorGroup.createAttribute("L",PredType::NATIVE_DOUBLE, H5S_SCALAR));

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

    DataSet timeVec_ds(integratorGroup.createDataSet("t_points",
                                                     PredType::NATIVE_FLOAT, space_t));
    DataSet temporalFreqVec_ds(integratorGroup.createDataSet("w_points",
                                                             PredType::NATIVE_FLOAT, space_t));

    DataSet spatialVec_ds(integratorGroup.createDataSet("s_points",
                                                        PredType::NATIVE_FLOAT, space_s));
    DataSet spatialFreqVec_ds(integratorGroup.createDataSet("k_points",
                                                            PredType::NATIVE_FLOAT, space_s));

    timeVec_ds.write(timeVec.memptr(), PredType::NATIVE_FLOAT);
    temporalFreqVec_ds.write(temporalFreqVec.memptr(), PredType::NATIVE_FLOAT);
    spatialVec_ds.write(spatialVec.memptr(), PredType::NATIVE_FLOAT);
    spatialFreqVec_ds.write(spatialFreqVec.memptr(), PredType::NATIVE_FLOAT);
}



void OutputManager::writeStimulusProperties(const Stimulus* stimulus)
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
        vec2 k = gratingStimulus->kVec();
        double orientation = gratingStimulus->orientation(true);
        double spatialFreq = gratingStimulus->spatialFreq();
        double w = gratingStimulus->temporalFreq();

        Attribute C_a(stim.createAttribute("C",PredType::NATIVE_DOUBLE, H5S_SCALAR));
        Attribute mask_a(stim.createAttribute("mask", StrType(PredType::C_S1, 64), H5S_SCALAR));
        Attribute maskSize_a(stim.createAttribute("mask_size",PredType::NATIVE_DOUBLE, H5S_SCALAR));
        Attribute k_a(stim.createAttribute("spatial_freq",PredType::NATIVE_DOUBLE, H5S_SCALAR));
        Attribute orientation_a(stim.createAttribute("orientation",PredType::NATIVE_DOUBLE, H5S_SCALAR));
        Attribute kx_a(stim.createAttribute("kx",PredType::NATIVE_DOUBLE, H5S_SCALAR));
        Attribute ky_a(stim.createAttribute("ky",PredType::NATIVE_DOUBLE, H5S_SCALAR));
        Attribute w_a(stim.createAttribute("temporal_freq",PredType::NATIVE_DOUBLE, H5S_SCALAR));



        C_a.write(PredType::NATIVE_DOUBLE, &C);
        mask_a.write( StrType(PredType::C_S1, 64), (&mask)->c_str());
        maskSize_a.write(PredType::NATIVE_DOUBLE, &maskSize);
        k_a.write(PredType::NATIVE_DOUBLE, &spatialFreq);
        orientation_a.write(PredType::NATIVE_DOUBLE, &orientation);
        kx_a.write(PredType::NATIVE_DOUBLE, &k(0));
        ky_a.write(PredType::NATIVE_DOUBLE, &k(1));
        w_a.write(PredType::NATIVE_DOUBLE, &w);
    }
}

void OutputManager::writeStimulus(const Stimulus* stimulus,
                                  const bool fourierTransform)
{

    herr_t status = H5Eset_auto1(NULL, NULL);
    status = H5Gget_objinfo (m_output->getId(), "/stimulus", 0, NULL);

    if (!status == 0){
        writeStimulusProperties(stimulus);
    }

    // Write stimuli
    Group stim = m_output->openGroup("/stimulus");

    fcube realStim = conv_to<fcube>::from(stimulus->spatioTemporal());
    writeDataSet(realStim, &stim, "spatio_temporal");

    if(fourierTransform){
        fcube complexStim = conv_to<fcube>::from(real(stimulus->fourierTransform()));
        writeDataSet(complexStim, &stim, "fourier_transform");
    }

}

void OutputManager::writeResponse(const Neuron& neuron,
                                  const bool fourierTransform)
{

    fcube response = conv_to<fcube>::from(neuron.response());

    string cellGroupName = neuron.type();
    herr_t status = H5Eset_auto1(NULL, NULL);
    status = H5Gget_objinfo (m_output->getId(), cellGroupName.c_str(), 0, NULL);

    if (!status == 0){
        Group cellGroup = m_output->createGroup(cellGroupName);
        string type = neuron.type();
        Attribute type_a(cellGroup.createAttribute("type",StrType(PredType::C_S1,64), H5S_SCALAR));
        type_a.write( StrType(PredType::C_S1, 64), (&type)->c_str());
    }

    //write response:
    Group res = m_output->createGroup(cellGroupName+"/response");

    writeDataSet(response, &res, "spatio_temporal");

    if(fourierTransform){
        fcube responseFT = conv_to<fcube>::from(real(neuron.responseFT()));
        writeDataSet(responseFT, &res, "fourier_transform");
    }

}

void OutputManager::writeImpulseResponse(const Neuron& neuron,
                                         const bool fourierTransform)
{

    fcube impulseResponse = conv_to<fcube>::from(neuron.impulseResponse());

    string cellGroupName = neuron.type();
    herr_t status = H5Eset_auto1(NULL, NULL);
    status = H5Gget_objinfo (m_output->getId(), cellGroupName.c_str(), 0, NULL);

    if (!status == 0){
        Group cellGroup = m_output->createGroup(cellGroupName);
        string type = neuron.type();
        Attribute type_a(cellGroup.createAttribute("type",StrType(PredType::C_S1,64), H5S_SCALAR));
        type_a.write( StrType(PredType::C_S1, 64), (&type)->c_str());
    }

    //write impulse response:
    Group impRes = m_output->createGroup(cellGroupName+"/impulse_response");
    writeDataSet(impulseResponse, &impRes, "spatio_temporal");

    if(fourierTransform){
        fcube impulseResponseFT =
                conv_to<fcube>::from(real(neuron.impulseResponseFourierTransform()));
        writeDataSet(impulseResponseFT, &impRes, "fourier_transform");
    }

}



void OutputManager::writeDataSet(const fcube data, Group* group, string name)
{

    hsize_t dim[3] = {data.n_slices, data.n_rows, data.n_cols};
    DataSpace space(3, dim);
    DataSet dataset(group->createDataSet(name, PredType::NATIVE_FLOAT, space));

    dataset.write(data.memptr(), PredType::NATIVE_FLOAT);

}

