import os
import numpy as np
from astropy.cosmology import FlatLambdaCDM
from collections import defaultdict

from ..core import CS
from ..core import utils as UT

cosmo = FlatLambdaCDM(H0=67.7, Om0=0.2592, Tcmb0=2.725)

class MergerTree:

    def __init__(self, sim):
        
        self.sim = sim
        self.base = sim.base
    
    def load_merger_file(self):

        filename = filename = os.path.join(self.base, "sink_mergers.txt")

        dtype = [
            ("Time", float),
            ("Redshift", float),

            ("ID1", int),
            ("ID2", int),

            ("Type1", int),
            ("Type2", int),

            ("M1", float),
            ("M2", float),

            ("P1x", float), ("P1y", float), ("P1z", float),
            ("P2x", float), ("P2y", float), ("P2z", float),

            ("V1x", float), ("V1y", float), ("V1z", float),
            ("V2x", float), ("V2y", float), ("V2z", float),
        ]

        data = np.genfromtxt(filename, names=True, dtype=dtype)

        return data
    
    def build_merger_tree(self, type=3):

        data = self.load_merger_file()

        tree = defaultdict(list)

        filename = os.path.join(self.base, "preliminary_info.txt")
        snaps, times, redshifts = np.genfromtxt(filename, usecols=(0, 1, 2), skip_header=1, unpack=True)

        for row in data:
            
            if row["Type1"] != type and row["Type2"] != type:
                continue

            descendant = row["ID1"]

            merger = {
                "time": row["Time"],
                "redshift": row["Redshift"],
                "CosmicTime": cosmo.age(row["Redshift"]).to("Myr").value,

                "primary": row["ID1"],
                "secondary": row["ID2"],

                "M1": row["M1"],
                "M2": row["M2"],
                "Mr": row["M1"] + row["M2"],
                "Mchirp": UT.chirp_mass(row["M1"], row["M2"]),
                "q": row["M2"] / row["M1"],
                "snap": int(UT.find_snaps_around_merger(row["Redshift"], snaps, redshifts)[0])
            }

            tree[descendant].append(merger)
        
        for ID in tree:
            tree[ID] = sorted(tree[ID], key=lambda x: x["CosmicTime"])
        
        if type == 3:
            self.sim.writer.write_pickle("merger_tree.pkl", tree)
        elif type == 2:
            self.sim.writer.write_pickle("popii_merger_tree.pkl", tree)
        
        return dict(tree)

