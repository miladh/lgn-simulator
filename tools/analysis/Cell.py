import numpy as np
import h5py
import argparse as ap

class Cell:
    """
    Class for neurons
    """

    def __init__(self, h5_group):
        self.group = h5_group
        for attr in h5_group.attrs.keys():
            setattr(self, attr, h5_group.attrs[attr])
            # print attr, h5_group.attrs[attr]

    def read_property(self, property=None, space=None):
        if hasattr(self, str(property)) and hasattr(getattr(self, property), str(space)):
            return 1

        if(property==None and space==None):
            for prop in self.group.keys():
                spaces = self.group[str(prop)].keys()
                data = {}
                for sp in spaces:
                    values = np.array(self.group[prop][sp])
                    data[str(sp)] = np.array(values)
                if not hasattr(self, prop): setattr(self, prop, ap.Namespace(**data))
        elif(property==None):
            for prop in self.group.keys():
                data = {}
                values = np.array(self.group[prop][space])
                data[str(space)] = np.array(values)
                if not hasattr(self, prop): setattr(self, prop, ap.Namespace(**data))
        elif(space==None):
            spaces = self.group[str(property)].keys()
            data = {}
            for sp in spaces:
                values = np.array(self.group[property][sp])
                data[str(sp)] = np.array(values)
            if not hasattr(self, property): setattr(self, property, ap.Namespace(**data))
        else:
            data = {}
            values = np.array(self.group[property][space])
            data[str(space)] = np.array(values)
            if not hasattr(self, property): setattr(self, property, ap.Namespace(**data))



    def zero_out_small_values(self):
       for attr in dir(self):
           if not attr.startswith('__') and not callable(getattr(self,attr)):
               property_dict = getattr(self, attr)

               for item in property_dict:
                #    low_values_indices = property_dict[item]< 1.0e-10
                #    print property_dict[item].shape, low_values_indices.shape, item
                   property_dict[item][property_dict[item]< 1.0e-10] = 0.0

               setattr(self, attr, property_dict)

    def normalize(self):
        for attr in dir(self):
            if not attr.startswith('__') and not callable(getattr(self,attr)):
                property_dict = getattr(self, attr)

                for item in property_dict:
                    values = property_dict[item]
                    for i in range(values.shape[0]):
                        if not abs(values[i,:,:]).max()==0:
                            values[i,:,:] /= abs(values[i,:,:]).max()

                setattr(self, attr, property_dict)
