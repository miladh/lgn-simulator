import numpy as np
import yaml
import os, os.path

from analysis.Simulation import Simulation

def simulation_extractor(sims, attr, value, return_as_list=True):
    """
    extracts simulations where attribute
    attr=value

    Parameters
    ----------
    sims : list
        list of Simulation objects
    attr: str
        name of the attribute

    Returns
    -------
    list
        list of Simulation objects where attribute
        attr=value

    """
    extracted_sims = []
    for sim in sims:
        p = sim.get_attribute(attr)
        if(p==value):
            extracted_sims.append(sim)

    if(not return_as_list and len(extracted_sims)==1):
        return extracted_sims[0]
    else:
        return extracted_sims


def extract_unique_simulation_attrs(sims, attr,  return_as_array=True):
    """
    extracts unique values for attribute attr
    from a list of simulations

    Parameters
    ----------
    sims : list
        list of Simulation objects
    attr: str
        name of the attribute

    Returns
    -------
    array
        array of unique attr values in the sims

    """

    attr_values = set()
    for sim in sims:
        attr_values.add(sim.get_attribute(attr))

    if( not return_as_array and len(attr_values)==1):
        return list(attr_values)[0]
    else:
        return np.array(sorted(list(attr_values)))


def get_simulations(data_path):
    """
    returns list of simulations

    Parameters
    ----------
    config_file : str
        data path

    Returns
    -------
    list
        list of Simulation objects.

    """
    import h5py
    sims=[]
    for root, dirs, files in os.walk(data_path):
        for dir in dirs:
            dir_name = os.path.join(root, dir)
            setting_file = os.path.join(dir_name, "_"+dir+".yaml")
            data_file = os.path.join(dir_name, dir+".h5")
            f = h5py.File(data_file, "r")
            sims.append(Simulation(setting_file, f))
    return sims
