import os
import numpy as np
from mpi4py import MPI
from collections import defaultdict

from ..physics.merger_tree import MergerTree
from ..core.constants import CS

class MergerSurroundings:
    
    def __init__(self, sim):
        
        self.sim = sim
        self.merger_tree = sim.load_mergertree()

        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()
        
        snaps = []
        primary = []
        secondary = []
        for ID in self.merger_tree:
            merger_branch = self.merger_tree[ID]
            for merger in merger_branch:
                primary.append(ID)
                secondary.append(merger["secondary"])
                snaps.append(merger["snap"])
        
        self.snaps = np.array(snaps)
        self.primary = np.array(primary)
        self.secondary = np.array(secondary)

    
    def find_surroundings(self):

        snaps = np.unique(self.snaps)
        surroundings = defaultdict(list)
        for i in range(self.rank, len(snaps), self.size):
            
            
            snap = snaps[i]
            print(f"Opening Simulation {self.sim.label}:Snapshot {snap}.")

            yt_utils = self.sim.load_yt(snap)
            ds = yt_utils.ds
            ad = ds.all_data()
            yt_utils.region_tools.add_metallicity()

            Galaxies = self.sim.reader.read_pickle(f"Galaxies_{snap:03d}.pkl")
            self.sim.binary_reader.read_sink_snap(snap)
            z = 1 / self.sim.binary_reader.time - 1
            Pos5 = self.sim.binary_reader.extract_data("Pos")
            ID5 = self.sim.binary_reader.extract_data("ID")
            Type = self.sim.binary_reader.extract_data("Type")

            ind = self.snaps == snap
            ID1 = self.primary[ind]
            ID2 = self.secondary[ind]

            for primary, secondary in zip(ID1, ID2):
                if primary not in ID5 or secondary not in ID5: continue
                data_already_present = False
                mask1 = primary == ID5
                p1 = Pos5[mask1][0]

                mask2 = secondary == ID5
                p2 = Pos5[mask2][0]

                dist = np.linalg.norm(p1-p2, axis=0) * 1e3 / (1 + z) / CS.hubble_parameter
                for i, Galaxy in Galaxies.items():
                    if primary == Galaxy["most_massive_id"] or primary in Galaxy["bh_ids"] or secondary == Galaxy["most_massive_id"] or secondary in Galaxy["bh_ids"]:

                        surrounding = {
                            "primary": primary,
                            "secondary": secondary,

                            "galaxy_pos": Galaxy["pos"],
                            "galaxy_mass": Galaxy["mass"],
                            "galaxy_radius": Galaxy["radius"],
                            "stellar_mass": Galaxy["stellar_mass"],
                            "gas_mass": Galaxy["gas_mass"],
                            "metallicity": Galaxy["metallicity"],
                            "separation": dist
                        }

                        surroundings[primary].append(surrounding)
                    
                        data_already_present = True
                
                
                if data_already_present == True:
                    continue

                region = yt_utils.region(p1, 0.1, units="kpc")
                mass2 = region[("PartType2", "Masses")].in_units("Msun")

                if np.sum(mass2) > 0:
                    continue
                
                bhmass, stellar, bhids, bhmasses, sep, mmax, mid = yt_utils.region_tools.FindBHandStellarProperty(ID5, Type)
                gasmass, coldgasmass, hotgasmass = yt_utils.region_tools.FindGasMass()

                Z = region[("PartType0", "Metallicity")]
                Zavg = np.mean(Z)

                surrounding = {
                    "primary": primary,
                    "secondary": secondary,

                    "galaxy_pos": None,
                    "galaxy_mass": None,
                    "galaxy_radius": None,
                    "stellar_mass": stellar,
                    "gas_mass": gasmass,
                    "metallicity": Zavg,
                    "separation": dist
                }

                surroundings[primary].append(surrounding)
                

        Alldicts = self.comm.gather(surroundings, root=0)

        if self.rank == 0:
            print("Gathered all surrounding data on rank 0.")
            Surroundings = defaultdict(list)

            for surroundings in Alldicts:
                for ID, surrounding in surroundings.items():
                    Surroundings[ID].extend(surrounding)
        
            self.sim.writer.write_pickle("merger_surroundings.pkl", dict(Surroundings))

            

            