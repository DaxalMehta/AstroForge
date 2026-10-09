import numpy as np
from dataclasses import dataclass, field, fields

from ..core.constants import CS

@dataclass
class Galaxy:

    # Indetity
    snapshot: int
    fof_index: int
    subhalo_index: int

    #Geometry
    position: np.ndarray | None = None
    velocity: np.ndarray | None = None
    aperture_radius: float | None = None
    center_of_mass: np.ndarray | None = None

    #Halo 
    catalog_mass: float | None = None
    catalog_mass_type: np.ndarray | None = None
    halfmass_radius: float | None = None
    halfmass_radius_type: np.ndarray | None = None
    vmax: float | None = None
    vmax_radius: float | None = None
    
    spin: float | None = None
    velocity_disperion: float | None = None

    #dark_matter_mass:
    dark_matter_mass: float | None = None

    #Gas
    gas_mass: float | None = None
    gas_cell_masses: np.ndarray = field(default_factory=lambda: np.empty(0))
    gas_cell_temperature: np.ndarray = field(default_factory=lambda: np.empty(0))
    gas_cell_density: np.ndarray = field(default_factory=lambda: np.empty(0))
    gas_cell_metallicity: np.ndarray = field(default_factory=lambda: np.empty(0))
    gas_cell_positions: np.ndarray = field(default_factory=lambda: np.empty((0, 3)))
    gas_cell_velocities: np.ndarray = field(default_factory=lambda: np.empty((0, 3)))

    #Stars
    stellar_mass: float | None = None
    stellar_particle_types: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=int))
    stellar_particle_ids: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=int))
    stellar_particle_masses: np.ndarray = field(default_factory=lambda: np.empty(0))
    stellar_particle_lifetimes: np.ndarray = field(default_factory=lambda: np.empty(0))
    stellar_particle_formation_times: np.ndarray = field(default_factory=lambda: np.empty(0))
    stellar_particle_positions: np.ndarray = field(default_factory=lambda: np.empty((0, 3)))
    stellar_particle_velocities: np.ndarray = field(default_factory=lambda: np.empty((0, 3)))

    #Black Holes
    bh_mass: float | None = None
    primary_bh_mass: float | None = None
    primary_bh_id: int | None = None
    bh_particle_ids: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=int))
    bh_particle_masses: np.ndarray = field(default_factory=lambda: np.empty(0))
    bh_separations: np.ndarray = field(default_factory=lambda: np.empty(0))
    bh_particle_positions: np.ndarray = field(default_factory=lambda: np.empty((0, 3)))
    bh_particle_velocities: np.ndarray = field(default_factory=lambda: np.empty((0, 3)))


    @property
    def baryonic_mass(self):
        return self.gas_mass + self.stellar_mass + self.bh_mass
    
    @property
    def gas_fraction(self):
        return self.gas_mass / self.baryonic_mass

    @property
    def baryonic_fraction(self):
        return self.baryonic_mass / self.halo_mass
    
    @property
    def stellar_fraction(self):
        return self.stellar_mass / self.baryonic_mass
    
    @property
    def cold_gas_mass(self):
        mask = (self.gas_cell_temperature < 1e4) & (self.gas_cell_number_densities > 1)
        return self.gas_cell_masses[mask].sum()

    @property
    def hot_gas_mass(self):
        return self.gas_mass - self.cold_gas_mass
    
    @property
    def cold_gas_fraction(self):
        return self.cold_gas_mass / self.gas_mass
    
    @property
    def hot_gas_fraction(self):
        return self.hot_gas_mass / self.gas_mass
    
    @property
    def mean_metallicity(self):
        return np.average(self.gas_cell_metallicity, weights=self.gas_cell_masses)

    @property
    def mean_temperature(self):
        return np.average(self.gas_cell_temperature, weights=self.gas_cell_masses)

    @property
    def gas_velocity_dispersion(self):

        v_com = np.average(self.gas_cell_velocities, axis=0, weight=self.gas_cell_masses)

        sigma = np.sqrt(np.sum(self.gas_cell_masses * np.sum((self.gas_cell_velocities - v_com)**2, axis=1))/ np.sum(self.gas_cell_masses))

        return sigma

    @property
    def gas_cell_number_densities(self):
        return self.gas_cell_density / 1.22 / CS.m_p_cgs

    @property
    def has_bh(self):
        return self.bh_mass > 0

    @property
    def n_bhs(self):
        return len(self.bh_particle_ids)
    
    @property
    def n_stars(self):
        return len(self.stellar_particle_ids)

    @property
    def n_popiii(self):
        return len(self.stellar_particle_ids[self.stellar_particle_types == CS.POPIII])

    @property
    def n_new_stars(self):
        return len(self.stellar_particle_ids[self.stellar_particle_types == CS.POPII_YOUNG])

    @property
    def n_old_stars(self):
        return len(self.stellar_particle_ids[self.stellar_particle_types == CS.POPII_OLD])
    
    @property
    def classification(self):
        ## Can be Star Forming, Quenched, 
        return "Galaxy"

    @property
    def stellar_half_mass_radius(self):
        return self.halfmass_radius_type[4]
        '''
        if self.n_stars == 0:
            return None

        r = np.linalg.norm(self.stellar_particle_positions - self.position, axis=1)
        order = np.argsort(r)
        r = r[order]
        m = self.stellar_particle_masses[order]

        cumulative_mass = np.cumsum(m)

        return np.interp(0.5 * cumulative_mass[-1], cumulative_mass, r)
        '''

@dataclass
class GalaxyCatalog:
    
    snapshot: int
    galaxies: list[Galaxy] = field(default_factory=list)

    _lookup: dict[int, Galaxy] = field(init=False, repr=False, default_factory=dict)

    def __post_init__(self):
        self._lookup = {
            galaxy.subhalo_index: galaxy
            for galaxy in self.galaxies
        }

    def __len__(self):
        return len(self.galaxies)

    def __iter__(self):
        return iter(self.galaxies)

    def __getitem__(self, item):
        return self.galaxies[item]

    def add(self, galaxy):
        self.galaxies.append(galaxy)

    @property
    def count(self):
        return len(self)

    @property
    def total_stellar_mass(self):
        return sum(g.stellar_mass for g in self)
    
    @property
    def stellar_masses(self):
        return np.array([g.stellar_mass for g in self])
    
    @property
    def largest_bh_masses(self):
        return np.array([g.primary_bh_mass for g in self])
    
    @property
    def largest_bh_ids(self):
        return np.array([g.primary_bh_id for g in self])

    @property
    def total_gas_mass(self):
        return sum(g.gas_mass for g in self)
    
    @property
    def total_bh_mass(self):
        return sum(g.bh_mass for g in self)
    
    @property
    def total_baryonic_mass(self):
        return sum(g.baryonic_mass for g in self)
    
    def with_black_holes(self):

        return GalaxyCatalog(
            snapshot=self.snapshot,
            galaxies=[g for g in self if g.has_bh]
        )

    def find_subhalo(self, subhalo_index: int) -> Galaxy:
        """
        Return the galaxy with the given subhalo index.

        Parameters
        ----------
        subhalo_index : int
            Subhalo index of the desired galaxy.

        Returns
        -------
        Galaxy
            The matching galaxy.

        Raises
        ------
        KeyError
            If no galaxy with the requested subhalo index exists.
        """

        #for galaxy in self.galaxies:
        #    if galaxy.subhalo_index == subhalo_index:
        #        return galaxy

        #raise KeyError(f"No galaxy with subhalo_index={subhalo_index}.")
        try:
            return self._lookup[subhalo_index]
        except KeyError:
            raise KeyError(f"No galaxy with subhalo_index={subhalo_index}.")









    def to_dict(self):
        """
        Convert the galaxy catalog to a dictionary of Numpy arrays
        suitable for serialization.
        """

        data = {"snapshot": self.snapshot}

        for field in fields(Galaxy):
            data[field.name] = [getattr(galaxy, field.name) for galaxy in self.galaxies]
        
        return data
    
    @classmethod
    def from_dict(cls, data):
        """
        Construct a GalaxyCatalog from a dictionary produced
        by 'to_dict()'
        """

        snapshot = data["snapshot"]

        field_names = {field.name for field in fields(Galaxy)}

        n_galaxies = len(data[next(iter(field_names))])

        galaxies = []
        for i in range(n_galaxies):
            galaxy_data = {
                name: data[name][i] for name in field_names
            }

            galaxies.append(Galaxy(**galaxy_data))
        
        return cls(snapshot=snapshot, galaxies=galaxies)