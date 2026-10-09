import numpy as np
import os
import yt
from dataclasses import dataclass, field, fields

from .particle import ParticleCatalog
from .galaxy import GalaxyCatalog
from ..inout.groups_reader import GroupReader
from ..visualizations.yt_utils import YTSession
from ..inout.binary_reader import BinaryReaderSnapshot

class Snapshot:
    
    def __init__(self, sim, snap):

        self.sim = sim
        self.snap = snap

        self._particles = None
        self._groups = None
        self._galaxies = None
        self._yt = None
        self._sink_snap = None
    
    @property
    def particles(self):
        
        if self._particles == None:
            self._particles = ParticleCatalog(self)
        
        return self._particles
    
    @property
    def groups(self):

        if self._groups == None:
            self._groups = GroupReader(self.sim.base, self.snap).read()
        return self._groups
    
    @property
    def galaxies(self):  
        
        if self._galaxies == None:
            data = self.sim.reader.read_pickle(f"Galaxies_{self.snap:03d}.pkl")
            self._galaxies = GalaxyCatalog.from_dict(data)
        return self._galaxies

    @property
    def yt(self):

        if self._yt == None:
            file = os.path.join(self.sim.base, f"snapdir_{self.snap:03d}/snap_{self.snap:03d}.0.hdf5")
            if not os.path.exists(file):
                file = os.path.join(self.sim.base, f"snap_{self.snap:03d}.hdf5")
            ds = yt.load(file)
            self._yt = YTSession(ds)

        return self._yt
    
    @property
    def sink_snap(self):

        if self._sink_snap == None:
            self._sink_snap = BinaryReaderSnapshot(self.sim).read_sink_snap(self.snap)
        
        return self._sink_snap



