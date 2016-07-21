#include "outputmanager.h"
#include <unistd.h>

#include "stimuli/grating/grating.h"
#include "stimuli/grating/circlemaskgrating.h"
#include "stimuli/grating/cscirclemaskgrating.h"
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
    fvec temporalFreqVec = conv_to<fvec>::from(FFTHelper::fftShift(integrator.temporalFreqVec()));
    fvec spatialFreqVec = conv_to<fvec>::from(FFTHelper::fftShift(integrator.spatialFreqVec()));


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


    //Grating attributes:-------------------------------------------------------------------------
    if (const  Grating * gratingStimulus = dynamic_cast<const Grating*>(stimulus) ) {
        double k = gratingStimulus->spatialFreq();
        double w = gratingStimulus->temporalFreq();
        double C = gratingStimulus->contrast();
        double phase = gratingStimulus->phase(true);
        double orientation = gratingStimulus->orientation(true);
        vec2 kVec = gratingStimulus->kVec();
        string mask = gratingStimulus->mask();


        Attribute k_a(stim.createAttribute("spatial_freq",PredType::NATIVE_DOUBLE, H5S_SCALAR));
        Attribute w_a(stim.createAttribute("temporal_freq",PredType::NATIVE_DOUBLE, H5S_SCALAR));
        Attribute C_a(stim.createAttribute("C",PredType::NATIVE_DOUBLE, H5S_SCALAR));
        Attribute phase_a(stim.createAttribute("phase",PredType::NATIVE_DOUBLE, H5S_SCALAR));
        Attribute orientation_a(stim.createAttribute("orientation",PredType::NATIVE_DOUBLE, H5S_SCALAR));
        Attribute kx_a(stim.createAttribute("kx",PredType::NATIVE_DOUBLE, H5S_SCALAR));
        Attribute ky_a(stim.createAttribute("ky",PredType::NATIVE_DOUBLE, H5S_SCALAR));
        Attribute mask_a(stim.createAttribute("mask", StrType(PredType::C_S1, 64), H5S_SCALAR));

        k_a.write(PredType::NATIVE_DOUBLE, &k);
        w_a.write(PredType::NATIVE_DOUBLE, &w);
        C_a.write(PredType::NATIVE_DOUBLE, &C);
        phase_a.write(PredType::NATIVE_DOUBLE, &phase);
        orientation_a.write(PredType::NATIVE_DOUBLE, &orientation);
        kx_a.write(PredType::NATIVE_DOUBLE, &kVec(0));
        ky_a.write(PredType::NATIVE_DOUBLE, &kVec(1));
        mask_a.write( StrType(PredType::C_S1, 64), (&mask)->c_str());


        //Circle mask///////////////////////////////////////////////////////////////////////////////////////
        if (const  CircleMaskGrating *circleMask = dynamic_cast<const CircleMaskGrating*>(gratingStimulus)){
            double maskSize = circleMask->maskSize();
            Attribute maskSize_a(stim.createAttribute("mask_size",PredType::NATIVE_DOUBLE, H5S_SCALAR));
            maskSize_a.write(PredType::NATIVE_DOUBLE, &maskSize);
        }


        //CSCircle mask------------------------------------------------------------------------------------
        if (const  CSCircleMaskGrating *csMask = dynamic_cast<const CSCircleMaskGrating*>(gratingStimulus)){
            double k_s = csMask->surroundSpatialFreq();
            double w_s = csMask->surroundTemporalFreq();
            double C_s = csMask->surroundContrast();
            double phase_s = csMask->surroundPhase(true);
            double orientation_s = csMask->surroundOrientation(true);
            double maskSize_s = csMask->surroundMaskSize();
            double maskSize = csMask->maskSize();
            vec2 kVec_s = gratingStimulus->kVec();


            Attribute k_s_a(stim.createAttribute("spatial_freq_sur",PredType::NATIVE_DOUBLE, H5S_SCALAR));
            Attribute w_s_a(stim.createAttribute("temporal_freq_sur",PredType::NATIVE_DOUBLE, H5S_SCALAR));
            Attribute C_s_a(stim.createAttribute("C_sur",PredType::NATIVE_DOUBLE, H5S_SCALAR));
            Attribute phase_s_a(stim.createAttribute("phase_sur",PredType::NATIVE_DOUBLE, H5S_SCALAR));
            Attribute orientation_s_a(stim.createAttribute("orientation_sur",PredType::NATIVE_DOUBLE, H5S_SCALAR));
            Attribute kx_s_a(stim.createAttribute("kx_sur",PredType::NATIVE_DOUBLE, H5S_SCALAR));
            Attribute ky_s_a(stim.createAttribute("ky_sur",PredType::NATIVE_DOUBLE, H5S_SCALAR));
            Attribute maskSize_s_a(stim.createAttribute("mask_size_sur",PredType::NATIVE_DOUBLE, H5S_SCALAR));
            Attribute maskSize_a(stim.createAttribute("mask_size",PredType::NATIVE_DOUBLE, H5S_SCALAR));


            k_s_a.write(PredType::NATIVE_DOUBLE, &k_s);
            w_s_a.write(PredType::NATIVE_DOUBLE, &w_s);
            C_s_a.write(PredType::NATIVE_DOUBLE, &C_s);
            phase_s_a.write(PredType::NATIVE_DOUBLE, &phase_s);
            orientation_s_a.write(PredType::NATIVE_DOUBLE, &orientation_s);
            kx_s_a.write(PredType::NATIVE_DOUBLE, &kVec_s(0));
            ky_s_a.write(PredType::NATIVE_DOUBLE, &kVec_s(1));
            maskSize_s_a.write(PredType::NATIVE_DOUBLE, &maskSize_s);
            maskSize_a.write(PredType::NATIVE_DOUBLE, &maskSize);

        }
    }

}

void OutputManager::writeStimulus(const Stimulus* stimulus)
{

    Group group = createGroupIfNotExist("stimulus");

    // Write stimuli
    fcube stim = conv_to<fcube>::from(stimulus->spatioTemporal());
    writeDataSet(stim, &group, "spatio_temporal");
}


void OutputManager::writeStimulusFourierTransform(const Stimulus* stimulus)
{
    cx_cube fftShiftedfourierTransform = stimulus->fourierTransform();
    //fftShift
    for(int i = 0; i < int(fftShiftedfourierTransform.n_slices); i++){
        fftShiftedfourierTransform.slice(i) =
                FFTHelper::fftShift(fftShiftedfourierTransform.slice(i));
    }
    cx_fcube stimFT = conv_to<cx_fcube>::from(fftShiftedfourierTransform);

    Group group = createGroupIfNotExist("stimulus/fourier_transform");
    writeDataSet(real(stimFT), &group, "real");
    writeDataSet(imag(stimFT), &group, "complex");
}


void OutputManager::writeResponse(const Neuron& neuron)
{
    fcube response = conv_to<fcube>::from(neuron.response());

    //write response:
    Group group = createGroupIfNotExist(neuron.type());
    if(!group.attrExists("type")){
        Attribute type(group.createAttribute("type",StrType(PredType::C_S1,64), H5S_SCALAR));
        type.write( StrType(PredType::C_S1, 64), (neuron.type()).c_str());
    }
    Group res = createGroupIfNotExist(neuron.type()+"/response");
    writeDataSet(response, &res, "spatio_temporal");


}

void OutputManager::writeImpulseResponse(const Neuron& neuron)
{

    fcube impulseResponse = conv_to<fcube>::from(neuron.impulseResponse());

    //write impulse response:
    Group group = createGroupIfNotExist(neuron.type());
    if(!group.attrExists("type")){
        Attribute type(group.createAttribute("type",StrType(PredType::C_S1,64), H5S_SCALAR));
        type.write( StrType(PredType::C_S1, 64), (neuron.type()).c_str());
    }

    Group impRes = createGroupIfNotExist(neuron.type()+"/impulse_response");
    writeDataSet(impulseResponse, &impRes, "spatio_temporal");

}


void OutputManager::writeResponseFourierTransform(const Neuron& neuron)
{
    cx_cube fftShiftedfourierTransform = neuron.responseFourierTransform();
    for(int i = 0; i < int(fftShiftedfourierTransform.n_slices); i++){
        fftShiftedfourierTransform.slice(i) =
                FFTHelper::fftShift(fftShiftedfourierTransform.slice(i));
    }
    cx_fcube responseFT = conv_to<cx_fcube>::from(fftShiftedfourierTransform);

    //write fourier transform of response:
    Group group = createGroupIfNotExist(neuron.type());
    if(!group.attrExists("type")){
        Attribute type(group.createAttribute("type",StrType(PredType::C_S1,64), H5S_SCALAR));
        type.write( StrType(PredType::C_S1, 64), (neuron.type()).c_str());
    }

    Group res = createGroupIfNotExist(neuron.type()+"/response");
    res.createGroup("fourier_transform");
    writeDataSet(real(responseFT), &res, "fourier_transform/real");
    writeDataSet(imag(responseFT), &res, "fourier_transform/complex");
}


void OutputManager::writeImpulseResponseFourierTransform(const Neuron& neuron)
{
    cx_cube fftShiftedfourierTransform = neuron.impulseResponseFourierTransform();
    for(int i = 0; i < int(fftShiftedfourierTransform.n_slices); i++){
        fftShiftedfourierTransform.slice(i) =
                FFTHelper::fftShift(fftShiftedfourierTransform.slice(i));
    }
    cx_fcube impulseResponseFT=conv_to<cx_fcube>::from(fftShiftedfourierTransform);

    //write fourier transform of impulse response:
    Group group = createGroupIfNotExist(neuron.type());
    if(!group.attrExists("type")){
        Attribute type(group.createAttribute("type",StrType(PredType::C_S1,64), H5S_SCALAR));
        type.write( StrType(PredType::C_S1, 64), (neuron.type()).c_str());
    }

    Group impRes = createGroupIfNotExist(neuron.type()+"/impulse_response");
    impRes.createGroup("fourier_transform");
    writeDataSet(real(impulseResponseFT), &impRes, "fourier_transform/real");
    writeDataSet(imag(impulseResponseFT), &impRes, "fourier_transform/complex");
}




Group OutputManager::createGroupIfNotExist(const string groupName)
{
    herr_t status = H5Eset_auto1(NULL, NULL);
    status = H5Gget_objinfo (m_output->getId(), groupName.c_str(), 0, NULL);

    Group group;
    if (!status == 0){
        group = m_output->createGroup(groupName);
    }else{
        group = m_output->openGroup(groupName);
    }

    return group;
}



void OutputManager::writeDataSet(const fcube data, Group* group, string name)
{

    hsize_t dim[3] = {data.n_slices, data.n_rows, data.n_cols};
    DataSpace space(3, dim);
    DataSet dataset(group->createDataSet(name, PredType::NATIVE_FLOAT, space));

    dataset.write(data.memptr(), PredType::NATIVE_FLOAT);

}


