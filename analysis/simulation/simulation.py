import os
import yt
import numpy as np
from dataclasses import dataclass, field, fields

from ..inout import Reader, Writer, GroupReader
from ..core.utils import file_range
from ..core.utils import CS
from ..physics import SinkParticles
from ..physics import MergerTree
from ..visualizations import YTSession
from ..simulation.registry import SIMULATIONS

from ..builders.binary_builder import ParticleBuilder, ParticleBuildState
from ..inout.binary_reader import BinaryReaderSnapshot

from ..builders.merger_builder import MergerBuilder

from ..objects.snapshot import Snapshot
from ..objects.merger import MergerCatalog
from ..objects.particle import ParticleHistoryDatabase, ParticleHistory
from ..objects.metadata import ParticleDatabase

class Simulation:
    def __init__(self, SIM, iso=False):

        print(f"Simulation: {SIM}")

        SIM = SIMULATIONS[SIM]
        self.SIM = SIM
        
        self.base = SIM["base"]
        self.base2 = SIM["base2"]
        self.feedback = SIM["feedback"]
        self.array_length = SIM["array_length"]

        self.color = SIM["color"]
        self.label = SIM["label"]
        self.linestyle = SIM["linestyle"]
        self.marker = SIM["marker"]
        
        root = SIM["label"]
        self.output_dir = f"{root}/"
        image_folder = f"Images/"
        
        os.makedirs(image_folder, exist_ok=True)
        os.makedirs(os.path.join(image_folder, self.output_dir), exist_ok=True)

        ## basic files
        self.start_snap, self.end_snap = file_range(self.base, mode=1)

        ## readers
        self.reader = Reader(self)
        self.writer = Writer(self)
        self.binary_reader_snapshot = BinaryReaderSnapshot(self)
        self.sink_tools = SinkParticles(self)

        self.merger_catalog = None

        self._particle_products = None
                
        self._history = History(self)
        
        
        self.build_state = self.binary_build_state()

        if self.build_state is None:
            self._build_particles()

        else:
            sink_base = f"{self.base}sink_particle_info/"
            files = set(sorted(os.path.join(sink_base, file) for file in os.listdir(sink_base)))
            processed = set(sorted(self.build_state.files_processed))

            #print(files)
            #print("***********************************************************************")
            #print(processed)

            new_files = files - processed

            if new_files:
                print("New files detected, running incremental files.")
                self._build_particles()
        
        

        ## Creating initial preliminary files
        if iso == False:
            filename = os.path.join(self.base, "preliminary_info.txt")
            if os.path.exists(filename):
                self.snaps, self.times, self.redshifts = np.genfromtxt(filename, usecols=(0, 1, 2), skip_header=1, unpack=True)
                self.snaps = self.snaps.astype(int)
                if np.max(self.snaps) != self.end_snap - 1:
                    self.snaps_time_redshifts(int(np.max(self.snaps)) + 1, self.end_snap, mode="a")
                    self.snaps, self.times, self.redshifts = np.genfromtxt(filename, usecols=(0, 1, 2), skip_header=1, unpack=True)
                    self.snaps = self.snaps.astype(int)
            else: 
                self.snaps_time_redshifts(self.start_snap, self.end_snap)
                self.snaps, self.times, self.redshifts = np.genfromtxt(filename, usecols=(0, 1, 2), skip_header=1, unpack=True)
                self.snaps = self.snaps.astype(int)

            

    def snapshot(self, snap):
        return Snapshot(self, snap)

    def _build_particles(self):
        if self._particle_products is None:
            self._particle_products = ParticleBuilder(self).build()
        
        return self._particle_products
    
    
    def metadata(self):
        if os.path.exists(os.path.join(self.base, "metadata.pkl")):
            data = self.reader.read_pickle("metadata.pkl")
            metadata = ParticleDatabase.from_dict(data)

        else:
            print("Crearing the metadata base and the history database")
            metadata = self._build_particles()[0]

        return metadata 

    def binary_build_state(self):
        filename = "particle_build_state.pkl"
        
        if os.path.exists(os.path.join(self.base, filename)):
            data =  self.reader.read_pickle(filename)
            return ParticleBuildState.from_dict(data)

        return None


    @property
    def history(self):
        return self._history

    def mergers(self):
        if self.merger_catalog is not None:
            return self.merger_catalog

        if os.path.exists(os.path.join(self.base, "merger_catalog.pkl")):
            data = self.reader.read_pickle("merger_catalog.pkl")
            self.merger_catalog = MergerCatalog.from_dict(data)
        else:
            catalog = MergerBuilder(self).build()

            self.writer.write_pickle("merger_catalog.pkl", catalog.to_dict())
            
            self.merger_catalog = catalog

        return self.merger_catalog

    def load_sink_exclusion_list(self):
        filename = os.path.join(self.base, "sink_exclusion_list.txt")
        if os.path.exists(filename):
            list = np.genfromtxt(filename, usecols=0, unpack=True)
            return list
        else:
            return np.array([])

    def snaps_time_redshifts(self, start_snap, end_snap, mode="w"):
        
        if self.start_snap is None: 
            return
        
        filename = f"{self.base}preliminary_info.txt"
        with open(filename, mode) as f:
        
            if mode == "w":
                f.write("Snap Time Redshift\n")
            
            for snap in range(start_snap, end_snap):
                sink_snap = self.binary_reader_snapshot.read_sink_snap(snap)

                time = sink_snap.time
                
                redshift = 1 / time - 1

                f.write(f"{snap} {time} {redshift}\n")
           
    def create_pickles(self):
        #self.sink_tools.load_full(force_create=True)
        self.sink_tools.load_reduced()
        #self.sink_tools.load_mergers(force_create=True)
        #self.sink_tools.load_popii_mergers(force_create=True) 

class History:
    def __init__(self, sim):
        
        self.sim = sim
        self.base = sim.base

        self._black_holes = None
        self._stellar_particles = None

    @property
    def black_holes(self):
       
        if self._black_holes is not None:
            return self._black_holes
        
        if os.path.exists(os.path.join(self.base, "black_holes.pkl")):
            data = self.sim.reader.read_pickle("black_holes.pkl")
            self._black_holes = {
                int(pid): ParticleHistory.from_dict(history)
                for pid, history in data.items()
            }
            
        
        else:
            history = self.sim._build_particles()[1]
            self._black_holes = history.black_holes
        
        return self._black_holes
    
    @property
    def stellar_particles(self):
        if self._stellar_particles is not None:
            return self._stellar_particles
        
        if os.path.exists(os.path.join(self.base, "stellar_particles.pkl")):
            data = self.sim.reader.read_pickle("stellar_particles.pkl")

            self._stellar_particles = {
                int(pid): ParticleHistory.from_dict(history)
                for pid, history in data.items()
            }
        
        else:
            history = self.sim._build_particles()[1]
            self._stellar_particles = history.stellar_particles
        
        return self._stellar_particles     
    
    
        
        
'''
 ### Old code

 class SinkEvolution:
    def __init__(self, sim):
        self.sim = sim

        self._sink_particles = None
        self._popii_particles = None
        self._mergermap = None
        self._bh_merger_tree = None
        self._popii_merger_tree = None
        self._merger_surroundings = None
    
    @property
    def sink_particles(self):
        
        if self._sink_particles == None:
            self._sink_particles = self.sim.sink_tools.load_full()
        
        return self._sink_particles
    
    @property
    def popii_particles(self):

        if self._popii_particles == None:
            self._popii_particles = self.sim.sink_tools.load_popii()
        
        return self._popii_particles
    
    @property
    def mergermap(self):
        
        if self._mergermap == None:
            self._mergermap = self.sim.sink_tools.load_mergers()

        return self._mergermap
    
    @property
    def bh_merger_tree(self):
        
        if self._bh_merger_tree == None:
            filename = "merger_tree.pkl" 
            if os.path.exists(os.path.join(self.base, filename)):
                self._bh_merger_tree = self.reader.read_pickle(filename)
            else:
                self._bh_merger_tree = MergerTree(self).build_merger_tree(type=3)
        
        return self._bh_merger_tree

    @property
    def popii_merger_tree(self):

        if self._popii_merger_tree == None:
            filename = "popii_merger_tree.pkl"
            if os.path.exists(os.path.join(self.base, filename)):
                self._popii_merger_tree = self.reader.read_pickle(filename)
            else:
                self._popii_merger_tree = MergerTree(self).build_merger_tree(type=2)
        
        return self._popii_merger_tree
    
    @property
    def mergersurroundings(self):
        
        if self._merger_surroundings == None:
            filename = "merger_surroundings.pkl"
            if os.path.exists(os.path.join(self.base, filename)):
                self._merger_surroundings = self.reader.read_pickle(filename)
                return self._merger_surroundings
            else:
                raise ImportError(f"No file existing.")           
 
 '''   
