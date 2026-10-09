import numpy as np
import os
import pickle
import struct
from ..core.utils import Type, SNeType

class Reader:
    def __init__(self, sim):
        self.sim = sim
        self.base = sim.base
    
    def read_pickle(self, filename):
        pickle_file = os.path.join(self.base, filename)
        
        print(f"Opening file {self.sim.label}:{filename}")
        with open(pickle_file, "rb") as f:
            data = pickle.load(f)
        
        return data
    
    def read_text(self, filename):
        file = os.path.join(self.base, filename)
        with open(file, "r") as f:
            header = f.readline().strip().split()
            data = f.readlines()
            dict_data = {h : [] for h in header}
            for i, line in enumerate(data):
                parts = line.split(" ")
                for j, part in enumerate(parts):
                    try:
                        dict_data[header[j]].append(float(part))
                    except ValueError:
                        dict_data[header[j]].append(part)
            
        return dict_data

class Writer:
    def __init__(self, sim):
        self.sim = sim
        self.base = sim.base
    
    def write_pickle(self, filename, dataset):
        pickle_file = os.path.join(self.base, filename)

        print(f"Writing file {self.sim.label}:{filename}")
        with open(pickle_file, "wb") as f:
            pickle.dump(dataset, f)

