import os
import yt
import numpy as np
from mpi4py import MPI

from ..visualizations.yt_utils import YTSession

class Snapshots:

    def __init__(self, sim):
        
        self.sim = sim
        
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()
    
    def find_volume(self):

        snaps = range(self.sim.start_snap, self.sim.end_snap)
        

        for i in range(self.rank, len(snaps), self.size):

            snap = snaps[i]
            snapfile = f"{self.sim.base}/snapdir_{snap:03d}/snap_{snap:03d}.0.hdf5"
            ds = yt.load(snapfile)
            ad = ds.all_data()

            ytutils = self.sim.load_yt(snap)

            print(f"[rank {self.rank}] Finding high res volume for snap {snap}")

            #V1, V2, V3 = ytutils.calculate_high_res_volume()

            #print(V1, V2, V3)
