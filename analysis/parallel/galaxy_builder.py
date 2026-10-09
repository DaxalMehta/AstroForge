import os
import yt
import numpy as np
from mpi4py import MPI

from ..builders.galaxy_builder import CatalogBuilder

class GalaxyBuilderMPI:

    def __init__(self, sim):
        
        self.sim = sim
        self.builder = CatalogBuilder(sim)
        
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()

    def build_all(self):

        snaps = range(self.sim.start_snap, self.sim.end_snap)

        for i in range(self.rank, len(snaps), self.size):

            snap = snaps[i]

            print(f"[rank {self.rank}] building galaxy catalog for snap {snap}")

            self.builder.build_catalog(snap)

        