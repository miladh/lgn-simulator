import numpy as np
import h5py

class Cell:
    """
    Class for neurons used in simulations

    """

    def __init__(self, cell_group):
        for attr in cell_group.attrs.keys():
            setattr(self, attr, cell_group.attrs[attr])
            # print attr, cell_group.attrs[attr]

        cell_properties = cell_group.keys()


        for property in cell_properties:
            spaces = cell_group[str(property)].keys()
            data = {}
            for space in spaces:
                values = np.array(cell_group[property][space])
                data[str(space)] = np.array(values)
                # print property, space, data[str(space)].shape
            setattr(self, property, data)

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
