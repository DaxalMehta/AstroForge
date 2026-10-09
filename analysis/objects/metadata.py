import numpy as np
import warnings
from dataclasses import dataclass, field, fields
from ..core.utils import find_cosmic_time

@dataclass
class ParticleDatabase:
    
    ids: np.ndarray
    initial_masses: np.ndarray
    particle_types: np.ndarray

    formation_times: np.ndarray
    stellar_lifetimes: np.ndarray

    status: np.ndarray | None = None
    final_masses: np.ndarray | None = None

    _sorter: np.ndarray = field(init=False, repr=False)

    def __post_init__(self):
        self._sorter = np.argsort(self.ids)
    
    def _lookup(self, ids, return_valid=False):
        """
        Return indices into the metadata arrays corresponding to the
        requested particle IDs.

        Parameters
        ----------
        ids : int or array-like
            Particle ID(s) to look up.
        return_valid : bool, optional
            If True, also return a boolean mask indicating which IDs
            were found.

        Returns
        -------
        int or ndarray
            Index (or indices) into the metadata arrays.
        bool or ndarray, optional
            Valid mask if `return_valid=True`.
        """
        self.ids = np.asarray(self.ids)

        scalar = np.isscalar(ids)
        ids = np.atleast_1d(ids)

        idx = np.searchsorted(self.ids, ids, sorter=self._sorter)

        lookup = np.full(idx.shape, -1, dtype=int)
        inside = idx < len(self.ids)

        if np.any(inside):
            lookup[inside] = self._sorter[idx[inside]]

        valid = np.zeros(ids.shape, dtype=bool)
        if np.any(inside):
            valid[inside] = self.ids[lookup[inside]] == ids[inside]

        if not np.all(valid):
            missing = ids[~valid]
            warnings.warn(f"{len(missing)} particle IDs were not found in the metadata.")

        if scalar:
            if return_valid:
                return lookup[0], valid[0]
            return lookup[0]

        if return_valid:
            return lookup, valid

        return lookup

    def initial_mass(self, ids):
        self.initial_masses = np.asarray(self.initial_masses)
        return self.initial_masses[self._lookup(ids)]

    def lifetime(self, ids):
        self.stellar_lifetimes = np.asarray(self.stellar_lifetimes)
        indices, valid = self._lookup(ids, return_valid=True)
        lifetimes = np.full(len(ids), np.nan)
        lifetimes[valid] = self.stellar_lifetimes[indices[valid]]
        return lifetimes

    def formation_time(self, ids):
        self.formation_times = np.asarray(self.formation_times)
        indices, valid = self._lookup(ids, return_valid=True)
        formation = np.full(len(ids), np.nan)
        formation[valid] = self.formation_times[indices[valid]]
        return formation

    def particle_type(self, ids):
        self.particle_types = np.asarray(self.particle_types)
        return self.types[self._lookup(ids)]

    def age(self, ids, at_time):
        """
        Return stellar ages (yr) at the supplied simulation time.
        """

        formation_times = self.formation_time(ids)

        return (
            find_cosmic_time(at_time)
            - find_cosmic_time(formation_times)
        )
    
    def to_dict(self):
        """
        Convert the metadata catalog to a dictionary of Numpy arrays
        suitable for serialization.
        """

        return {
            f.name: getattr(self, f.name) 
            for f in fields(self) if f.init
        }
    
    @classmethod
    def from_dict(cls, data):
        """
        Construct a ParticleDatabase from a dictionary produced
        by `to_dict()`.
        """

        return cls(**data)