import numpy as np
from dataclasses import dataclass, field, fields
from scipy.spatial import cKDTree

from ..core.constants import CS
from ..core.enums import SinkType, ParticleType
from ..inout.reader_writer import Reader, Writer
from ..objects.metadata import ParticleDatabase


@dataclass
class BaseParticleSet:

    positions: np.ndarray
    velocities: np.ndarray
    masses: np.ndarray

    ids: np.ndarray | None = None
    
    tree: cKDTree = field(init=False, repr=False)

    def __post_init__(self):
        self.tree = cKDTree(self.positions)

    def _subset(self, mask):

        kwargs = {}

        n = self.count

        for f in fields(self):

            # skip the KDTree
            if f.name == "tree":
                continue

            value = getattr(self, f.name)

            #subset numpy arrays
            if (isinstance(value, np.ndarray) and len(value) == n):
                kwargs[f.name] = value[mask]
            
            # Leave everything else unchanged
            else:
                kwargs[f.name] = value
        
        return self.__class__(**kwargs)

    def query(self, center, radius):

        return self.tree.query_ball_point(center, radius)

    def inside(self, center, radius):

        idx = self.tree.query_ball_point(center, radius)
        return self._subset(idx)
    
    @property
    def total_mass(self):
        return self.masses.sum()
    
    @property
    def count(self):
        return len(self.masses)
    
    @property
    def center_of_mass(self):
        if self.count == 0:
            return None
        return np.average(self.positions, axis=0, weights=self.masses)
    
    @property
    def largest_particle(self):
        if self.count == 0:
            return None, 0.0
        i = np.argmax(self.masses)
        pid = None if self.ids is None else self.ids[i]
        return pid, self.masses[i]
    
    @property
    def most_massive(self):
        if self.count == 0:
            return None
        return self._subset([np.argmax(self.masses)])

@dataclass
class GasParticles(BaseParticleSet):
    sizes: np.ndarray | None = None
    densities: np.ndarray | None = None
    number_densities: np.ndarray | None = None
    temperatures: np.ndarray | None = None
    metallicities: np.ndarray | None = None

@dataclass
class ContaminationParticles(BaseParticleSet):
    pass

@dataclass
class StarParticles(BaseParticleSet):

    formation_times: np.ndarray | None = None
    stellar_lifetimes: np.ndarray | None = None

@dataclass
class Sinks(BaseParticleSet):

    particle_types: np.ndarray | None = None
    formation_times: np.ndarray | None = None
    stellar_lifetimes: np.ndarray | None = None
    
    @property
    def black_holes(self):
        return self._subset(self.particle_types == 3)
    
    @property
    def popIII(self):
        return self._subset(self.particle_types == 0)
    
    @property
    def popII(self):
        return self._subset(self.particle_types == 2)
    
    def distances_to(self, pos):
        return np.linalg.norm(self.positions - pos, axis=1)

class ParticleCatalog:
    def __init__(self, snapshot):

        self.sim = snapshot.sim
        self.snap = snapshot.snap
        
        self.binary_reader_snapshot = snapshot.sim.binary_reader_snapshot

        yt = snapshot.yt
        yt.fields.metallicity()

        self.ds = yt.ds
        self.ad = yt.ad

        self.contamination = self._load_contamination()
        self.stars = self._load_stars()
        self.sinks = self._load_sinks()
        self.gas = self._load_gas()
    
    def _load_contamination(self):       
        positions = self.ad[("PartType2", "Coordinates")].v
        velocities = self.ad[("PartType2", "Velocities")].in_units("km/s").v
        masses = self.ad[("PartType2", "Masses")].in_units("Msun").v

        return ContaminationParticles(positions=positions, velocities=velocities, masses=masses)
    
    def _load_gas(self):
    
        positions = self.ad[("PartType0", "Coordinates")].v
        velocities = self.ad[("PartType0", "Velocities")].in_units("km/s").v
        masses = self.ad[("PartType0", "Masses")].in_units("Msun").v

        densities = self.ad[("gas", "density")].in_units("g/cm**3").v
        number_densities = self.ad[("gas", "number_density")].v
        temperatures = self.ad[("gas", "temperature")].v
        metallicities = self.ad[("PartType0", "Metallicity")].v

        volumes = masses * CS.Msun_to_g / densities
        sizes = (3 * volumes / 4 / np.pi) ** (1/3)
        sizes /= CS.pc_cgs

        return GasParticles(positions=positions, velocities=velocities, masses=masses, sizes=sizes, densities=densities, number_densities=number_densities, temperatures=temperatures, metallicities=metallicities)
    
    
    def _load_stars(self):

        if "PartType4" not in self.ds.fields:
            return StarParticles(positions=np.empty((0, 3)), velocities=np.empty((0, 3)), masses=np.empty(0), ids=np.empty(0, dtype=np.int64), formation_times=np.empty(0), stellar_lifetimes=np.empty(0))

        positions = self.ad[("PartType4", "Coordinates")].v
        velocities = self.ad[("PartType4", "Velocities")].in_units("km/s").v
        masses = self.ad[("PartType4", "Masses")].in_units("Msun").v
        ids = self.ad[("PartType4", "ParticleIDs")].v.astype(int)

        data = self.sim.reader.read_pickle("metadata.pkl")
        metadata = ParticleDatabase.from_dict(data)
        #print(ids)
        formation_times = metadata.formation_time(ids)
        #print(formation_times)
        stellar_lifetimes = metadata.lifetime(ids)

        return StarParticles(positions=positions, velocities=velocities, masses=masses, ids=ids, formation_times=formation_times, stellar_lifetimes=stellar_lifetimes)

    def _load_sinks(self):

        if "PartType5" not in self.ds.fields:
            return Sinks(positions=np.empty((0, 3)), velocities=np.empty((0, 3)), masses=np.empty(0), ids=np.empty(0, dtype=np.int64), particle_types=np.empty(0, dtype=np.int64), formation_times=np.empty(0), stellar_lifetimes=np.empty(0))
        sink_snap = self.binary_reader_snapshot.read_sink_snap(self.snap)
        
        positions = sink_snap.positions
        velocities = sink_snap.velocities * pow(sink_snap.time, 0.5)
        masses = sink_snap.masses * 1e10 / CS.hubble_parameter
        ids = sink_snap.ids
        particle_types = sink_snap.particle_types
        formation_times = sink_snap.formation_times
        stellar_lifetimes = sink_snap.stellar_lifetimes

        return Sinks(positions=positions, velocities=velocities, masses=masses, ids=ids, particle_types=particle_types, formation_times=formation_times, stellar_lifetimes=stellar_lifetimes)

@dataclass
class ParticleHistory:
    """
    Complete time evolution of a single sink particle.
    """
    id: int
    particle_type: int

    dump_times: list[float] = field(default_factory=list)
    masses: list[float] = field(default_factory=list)
    merger_masses: list[float] = field(default_factory=list)

    positions: list[np.ndarray] = field(default_factory=list)
    velocities: list[np.ndarray] = field(default_factory=list)

    def append(self, record):
        """Append a new SinkInfoRecord to the history."""

        self.dump_times.append(record.dump_time)
        self.masses.append(record.mass)
        self.merger_masses.append(record.merger_mass)

        self.positions.append(record.position)
        self.velocities.append(record.velocity)
    
    @property
    def formation_time(self):
        return self.dump_times[0]
    
    @property
    def final_mass(self):
        return self.masses[-1]
    
    @property
    def initial_mass(self):
        return self.masses[0]
    
    @property
    def lifetime(self):
        return self.dump_times[-1] - self.dump_times[0]
    
    @property
    def n_outputs(self):
        return len(self.dump_times)
    
    def as_arrays(self):
        """
        Convert all internally stored lists to numpy arrays
        """

        self.dump_times = np.asarray(self.dump_times)
        self.masses = np.asarray(self.masses)
        self.merger_masses = np.asarray(self.merger_masses)
        
        self.positions = np.asarray(self.positions)
        self.velocities = np.asarray(self.velocities)
    
    def to_dict(self):
        """
        Convert the particle history to a dictionary suitable
        for serialization.
        """

        return {
            f.name: getattr(self, f.name) 
            for f in fields(self)
        }


    @classmethod
    def from_dict(cls, data):
        """
        Construct a ParticleHistory from a dictionary produced
        by `to_dict()`.
        """

        return cls(**data)

@dataclass        
class ParticleHistoryDatabase:
    
    """
    Collection of particle histories.
    """

    black_holes: dict[int, ParticleHistory] = field(default_factory=dict)
    stellar_particles: dict[int, ParticleHistory] = field(default_factory=dict)

    @property
    def n_black_holes(self):
        return len(self.black_holes)

    @property
    def n_stars(self):
        return len(self.stellar_particles)

    @property
    def n_particles(self):
        return self.n_black_holes + self.n_stars

    @property
    def black_hole_ids(self):
        return list(self.black_holes.keys())

    @property
    def stellar_ids(self):
        return list(self.stellar_particles.keys())

    def as_arrays(self):
        """
        Convert all particle histories to NumPy arrays.
        """

        for history in self.black_holes.values():
            history.as_arrays()

        for history in self.stellar_particles.values():
            history.as_arrays()

