import numpy as np
import os
from collections import defaultdict
from dataclasses import dataclass, field, fields
from ..core.utils import chirp_mass, find_cosmic_time, find_snaps_around_merger

@dataclass(slots=True)
class Merger:
    time: float
    redshift: float

    primary_id: int
    secondary_id: int

    primary_type: int
    secondary_type: int

    primary_mass: float
    secondary_mass: float

    snapshot: int

    snr_lisa: float | None = None
    snr_et: float | None = None
    snr_lgwa: float | None = None

    @property
    def remnant_mass(self):
        return self.primary_mass + self.secondary_mass
    
    @property
    def chirp_mass(self):
        return chirp_mass(self.primary_mass, self.secondary_mass)
    
    @property
    def mass_ratio(self):
        return self.secondary_mass / self.primary_mass
    
    @property
    def cosmic_time(self):
        return find_cosmic_time(self.time)

@dataclass(slots=True)    
class MergerCatalog:
    
    mergers: list[Merger]

    def __len__(self):
        return len(self.mergers)
    
    def __iter__(self):
        return iter(self.mergers)
    
    def __getitem__(self, index):
        return self.mergers[index]

    def __repr__(self):
        return f"MergerCatalog({len(self)} mergers)"

    @property
    def primary_ids(self):
        return np.array([m.primary_id for m in self.mergers])
    
    @property
    def primary_masses(self):
        return np.array([m.primary_mass for m in self.mergers])
    
    @property
    def secondary_masses(self):
        return np.array([m.secondary_mass for m in self.mergers])

    @property
    def binary_masses(self):
        return np.array([m.primary_mass for m in self.mergers]), np.array([m.secondary_mass for m in self.mergers])
    
    @property
    def snapshots(self):
        return np.array([m.snapshot for m in self.mergers])
    
    @property
    def get_merger_arrays(self):

        m1 = np.array([m.primary_mass for m in self.mergers])
        m2 = np.array([m.secondary_mass for m in self.mergers])
        z = np.array([m.redshift for m in self.mergers])
        snr_lisa = np.array([m.snr_lisa for m in self.mergers])
        snr_et = np.array([m.snr_et for m in self.mergers])
        snr_lgwa = np.array([m.snr_lgwa for m in self.mergers])

        return m1, m2, z, snr_lisa, snr_et, snr_lgwa
    
    def subset(self, *, primary_type=None, secondary_type=None, primary_id=None, secondary_id=None, snapshot=None, redshift_range=None, primary_mass_range=None, secondary_mass_range=None):
        """
        Return a subset of the merger catalog satisfying the given criteria.

        Parameters
        ----------
        primary_type : int or SinkType, optional
        secondary_type : int or SinkType, optional
        primary_id : int, optional
        secondary_id : int, optional
        snapshot : int, optional
        redshift_range : tuple(float, float), optional
            (z_min, z_max)
        primary_mass_range : tuple(float, float), optional
            (m_min, m_max)
        secondary_mass_range : tuple(float, float), optional
            (m_min, m_max)
        """

        mergers = self.mergers

        if primary_type is not None:
            mergers = [m for m in mergers if m.primary_type == primary_type]

        if secondary_type is not None:
            mergers = [m for m in mergers if m.secondary_type == secondary_type]

        if primary_id is not None:
            mergers = [m for m in mergers if m.primary_id == primary_id]

        if secondary_id is not None:
            mergers = [m for m in mergers if m.secondary_id == secondary_id]

        if snapshot is not None:
            mergers = [m for m in mergers if m.snapshot == snapshot]

        if redshift_range is not None:
            zmin, zmax = redshift_range
            mergers = [m for m in mergers if zmin <= m.redshift <= zmax]

        if primary_mass_range is not None:
            mmin, mmax = primary_mass_range
            mergers = [m for m in mergers if mmin <= m.primary_mass <= mmax]

        if secondary_mass_range is not None:
            mmin, mmax = secondary_mass_range
            mergers = [m for m in mergers if mmin <= m.secondary_mass <= mmax]

        return MergerCatalog(mergers)

    def merger_lookup(self, type=3):
        
        tree = defaultdict(list)

        for merger in self.mergers:
            if merger.primary_type != type:
                continue
            tree[merger.primary_id].append(merger)
        
        for ID in tree:
            tree[ID] = sorted(tree[ID], key=lambda x: x.time)

        return dict(tree)


    def to_dict(self):
        """
        Convert the merger catalog to a dictionary of Numpy arrays 
        suitable for saving to a pickle file.
        """    

        return {
            field.name: np.array([getattr(m, field.name) for m in self.mergers])
            for field in fields(Merger)
        }

    @classmethod
    def from_dict(cls, data):
        """
        Construct a MergerCatalog from a dictionary produced 
        by 'to_dict()'.
        """

        mergers = [
            Merger(**{name: data[name][i] for name in data})
            for i in range(len(next(iter(data.values()))))
        ]

        return cls(mergers)
