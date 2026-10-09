import numpy as np
from dataclasses import dataclass, field, fields

"""
<KeysViewHDF5 ['GroupFirstSub', 'GroupLen', 'GroupLenType', 'GroupMass', 'GroupMassType', 'GroupNsubs', 'GroupOffsetType', 'GroupPos', 'GroupVel', 'Group_M_Crit200', 'Group_M_Crit500', 'Group_M_Mean200', 'Group_M_TopHat200', 'Group_R_Crit200', 'Group_R_Crit500', 'Group_R_Mean200', 'Group_R_TopHat200']> 

<KeysViewHDF5 ['SubhaloCM', 'SubhaloGroupNr', 'SubhaloHalfmassRad', 'SubhaloHalfmassRadType', 'SubhaloIDMostbound', 'SubhaloIntertiaTensorStars', 'SubhaloKineticEnergyStars', 'SubhaloLen', 'SubhaloLenType', 'SubhaloMass', 'SubhaloMassInHalfRad', 'SubhaloMassInHalfRadType', 'SubhaloMassInMaxRad', 'SubhaloMassInMaxRadType', 'SubhaloMassInRad', 'SubhaloMassInRadType', 'SubhaloMassType', 'SubhaloOffsetType', 'SubhaloParentRank', 'SubhaloPos', 'SubhaloRankInGr', 'SubhaloRotationalEnergyStars', 'SubhaloSpin', 'SubhaloSpinType', 'SubhaloVel', 'SubhaloVelDisp', 'SubhaloVmax', 'SubhaloVmaxRad']>
"""

HALO_FIELDS = {
    "GroupFirstSub": "first_subhalo",
    "GroupLen": "n_particles",
    "GroupLenType": "n_particles_type",
    "GroupMass": "mass",
    "GroupMassType": "mass_type",
    "GroupNsubs": "n_subs",
    "GroupOffsetType": "offset_type",
    "GroupPos": "position",
    "GroupVel": "velocity",
    "Group_M_Crit200": "mass_200",
    "Group_M_Crit500": "mass_500",
    "Group_M_Mean200": "mass_mean_200",
    "Group_M_TopHat200": "mass_top_200",
    "Group_R_Crit200": "radius_200",
    "Group_R_Crit500": "radius_500",
    "Group_R_Mean200": "radius_mean_200",
    "Group_R_TopHat200": "radius_top_200",
}

SUBHALO_FIELDS = {
    "SubhaloCM": "center_of_mass",
    "SubhaloGroupNr": "group_index",
    "SubhaloHalfmassRad": "halfmass_radius",
    "SubhaloHalfmassRadType": "halfmass_radius_type",
    "SubhaloIDMostbound": "most_bound_id",
    "SubhaloIntertiaTensorStars": "inertia_tensor_stars",
    "SubhaloKineticEnergyStars": "kinetic_energy_stars",
    "SubhaloLen": "n_particles",
    "SubhaloLenType": "n_particles_type",
    "SubhaloMass": "mass",
    "SubhaloMassInHalfRad": "mass_in_half_radius",
    "SubhaloMassInHalfRadType": "mass_in_half_radius_type",
    "SubhaloMassInMaxRad": "mass_in_max_radius",
    "SubhaloMassInMaxRadType": "mass_in_max_radius_type",
    "SubhaloMassInRad": "mass_in_radius",
    "SubhaloMassInRadType": "mass_in_radius_type",
    "SubhaloMassType": "mass_type",

    "SubhaloOffsetType": "offset_type",
    "SubhaloParentRank": "parent_subhalo_rank",
    "SubhaloPos": "position",
    "SubhaloRankInGr": "rank_in_group",
    "SubhaloRotationalEnergyStars": "rotation_energy_stars",
    "SubhaloSpin": "spin",
    "SubhaloSpinType": "spin_type",
    "SubhaloVel": "velocity",
    "SubhaloVelDisp": "velocity_dispersion",
    "SubhaloVmax": "vmax",
    "SubhaloVmaxRad": "vmax_radius",
}
@dataclass
class Header:
    pass

    @property
    def redshift(self):
        return self["Redshift"]

@dataclass
class HaloCatalog:
    pass

    @property
    def stellar_mass(self):
        return self.mass_type[:, 4]
    
    @property
    def gas_mass(self):
        return self.mass_type[:, 0]

@dataclass
class SubhaloCatalog:
    pass

    @property
    def stellar_mass(self):
        return self.mass_type[:, 4]

    @property
    def gas_mass(self):
        return self.mass_type[:, 0]

@dataclass
class Groups:

    snapshot: int
    header: Header
    halos: HaloCatalog
    subhalos: SubhaloCatalog





