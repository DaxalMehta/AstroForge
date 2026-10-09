from dataclasses import dataclass
import numpy as np

from ..core.enums import SinkType

@dataclass(slots=True)
class SinkInfoRecord:
    """
    One record from a sink_particle_info file.
    """
    dump_time: float

    position: np.ndarray
    velocity: np.ndarray
    acceleration: np.ndarray

    mass: float
    formation_time: float
    
    id: int

    particle_type: int

    stellar_lifetime: float
    merger_mass: float

@dataclass(slots=True)
class SinkDump:
    """
    One output (dump) from a sink_particle_info file.
    """

    time: float
    records: list[SinkInfoRecord]

@dataclass(slots=True)
class SinkSnapshotRecord:
    """
    One record from a sink_snap file.
    """

    position: np.ndarray
    velocity: np.ndarray
    acceleration: np.ndarray

    mass: float
    formation_time: float

    id: int

    particle_type: int

    stellar_lifetime: float

    merger_mass: float

    number_of_sne: int | None = None

@dataclass(slots=True)
class SinkSnapshot:
    """
    The entire sink data processed from sink_snap
    """

    time: float
    redshift: float

    sink_particles: list[SinkSnapshotRecord]

    @property
    def num_sinks(self):
        return len(self.sink_particles)

    @property
    def masses(self):
        return np.array([sink_particle.mass for sink_particle in self.sink_particles])

    @property
    def ids(self):
        return np.array([sink_particle.id for sink_particle in self.sink_particles])
    
    @property
    def particle_types(self):
        return np.array([sink_particle.particle_type for sink_particle in self.sink_particles])
    
    @property
    def positions(self):
        return np.array([sink_particle.position for sink_particle in self.sink_particles])
    
    @property
    def velocities(self):
        return np.array([sink_particle.velocity for sink_particle in self.sink_particles])
    
    @property
    def num_snes(self):
        return np.array([sink_particle.number_of_sne for sink_particle in self.sink_particles])
    
    @property
    def formation_times(self):
        return np.array([sink_particle.formation_time for sink_particle in self.sink_particles])
    
    @property
    def stellar_lifetimes(self):
        return np.array([sink_particle.stellar_lifetime for sink_particle in self.sink_particles])