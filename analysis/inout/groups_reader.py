import illustris_python as il
import numpy as np
import os
import h5py
import glob
from hmf import MassFunction
from astropy.cosmology import Planck15 as cosmo
from collections import defaultdict

from ..core import CS
from ..objects.groups import (Groups, HaloCatalog, SubhaloCatalog, Header, HALO_FIELDS, SUBHALO_FIELDS)


class GroupReader:
    """
    Read FoF/Subfind catalogs from AREPO.

    Parameters
    ----------
    base : str
        Simulation directory.
    snap : int
        Snapshot number.
    """

    def __init__(self, base, snap):
        self.base = base
        self.snap = snap

        self.files = self._find_files()

    def read(self):

        header = self._read_header()
        halos = self._read_catalog("Group", HALO_FIELDS, HaloCatalog)
        subhalos = self._read_catalog("Subhalo", SUBHALO_FIELDS, SubhaloCatalog)

        return Groups(snapshot=self.snap,
                      header=header,
                      halos=halos,
                      subhalos=subhalos)

    def _find_files(self):
        """
        Return all FoF/Subfind files belonging to this snapshot.
        """

        directory = os.path.join(self.base, f"groups_{self.snap:03d}")

        files = sorted(glob.glob(os.path.join(directory, f"fof_subhalo_tab_{self.snap:03d}.*.hdf5")))

        if len(files) == 0:
            raise FileNotFoundError(f"No FoF catalog found for snapshot {self.snap}.")

        return files
    
    def _read_header(self):
        header = Header()

        with h5py.File(self.files[0], "r") as f:
            for key, value in f["Header"].attrs.items():
                setattr(header, key, value)

        return header

    def _read_catalog(self, group_name, field_map, cls):
        """
        Read one catalog (Group or Subhalo).
        """

        arrays = {}

        for filename in self.files:
            with h5py.File(filename, "r") as f:

                if group_name not in f: 
                    continue
                group = f[group_name]

                for hdf_name, attr_name in field_map.items():
                    if hdf_name not in group:
                        continue
                    arrays.setdefault(attr_name, []).append(group[hdf_name][:])

        catalog = cls()

        for name, chunks in arrays.items():
            setattr(catalog, name, np.concatenate(chunks, axis=0))

        return catalog
        

'''
class GroupReader:
    def __init__(self, base, snap, halo_fields=None, subhalo_fields=None):
        self.base = base
        self.snap = snap
        self.halo_fields = ["GroupMass", "GroupPos", "Group_R_TopHat200", "GroupMassType"]
        if halo_fields:
            self.halo_fields.extend(halo_fields)
        self.subhalo_fields = ['SubhaloMass','SubhaloCM', 'SubhaloHalfmassRad', 'SubhaloGroupNr', 'SubhaloMassType', 'SubhaloOffsetType', 'SubhaloLenType']
        if subhalo_fields:
            self.subhalo_fields.extend(subhalo_fields)
        
        self.halos = self._ReadHalos()
        self.subhalos = self._ReadSubHalos()
        self.header = self._ReadHeader()
        self.Redshift = self.header['Redshift']
        
        if "GroupPos" not in self.halos:
            return
        self.HaloPos = self.halos["GroupPos"]
        self.HaloRad = self.halos["Group_R_TopHat200"]
        self.HaloMassType = self.halos["GroupMassType"]
        if "SubhaloCM" not in self.subhalos:
            return
        self.SubhaloPos = self.subhalos["SubhaloCM"]
        self.SubhaloRad = self.subhalos["SubhaloHalfmassRad"]
        self.SubhaloMass = self.subhalos["SubhaloMass"]
        self.SubhaloMassType = self.subhalos["SubhaloMassType"]
        self.SubhaloGrNr = self.subhalos["SubhaloGroupNr"]
        self.SubhaloOffsetType = self.subhalos["SubhaloOffsetType"]
        self.SubhaloLenType = self.subhalos["SubhaloLenType"]

    def _ReadHeader(self):
        header = il.groupcat.loadHeader(self.base, self.snap)
        return header
    
    def _ReadHalos(self):
        halos = il.groupcat.loadHalos(self.base, self.snap, self.halo_fields)
        return halos
    
    def _ReadSubHalos(self):
        subhalos = il.groupcat.loadSubhalos(self.base, self.snap, self.subhalo_fields)
        return subhalos
    
    @property
    def boxsize(self):
        return self.header["BoxSize"]

    def HostHaloProperties(self, SinkPos):
        dist = np.linalg.norm(self.HaloPos - SinkPos, axis=1)
        ind = np.argmin(dist)
        if (self.HaloMassType[ind, 5] == 0): 
            return None, None
        else:
            return self.HaloPos[ind], self.HaloRad[ind]
    
    def HostSubhaloProperties(self, SinkPos):
        dist = np.linalg.norm(self.SubhaloPos - SinkPos, axis=1)
        ind = np.argmin(dist)
        #if (self.SubhaloMassType[ind, 5] == 0):
        #    return None, None
        #else:
        return self.SubhaloPos[ind], self.SubhaloRad[ind]
    
    def HostGalaxyProperties(self, SinkPos):
        
        SubhaloMass = self.subhalos["SubhaloMass"] * 1e10 / CS.hubble_parameter
        SubhaloPos = self.SubhaloPos
        SubhaloRad = self.SubhaloRad
        
        dist = np.linalg.norm(SubhaloPos - SinkPos, axis=1)
        mask = dist < SubhaloRad
        
        if not np.any(mask):
            return None, None, None
        
        valid_inds = np.where(mask)[0]
        best = valid_inds[np.argmax(SubhaloMass[valid_inds])]
        
        return best, SubhaloPos[best], SubhaloRad[best]
    
    def Normalization(self, Volume):
        HaloMasses = self.SubhaloMass * 1e10 
        SinkMasses = self.SubhaloMassType[:, 5] * 1e10 
        
        mask = SinkMasses > 0
        if mask.sum() == 0:
            return None
        
        mass_low = np.min(HaloMasses[mask])
        mass_high = np.max(HaloMasses[mask])

        print(f"Low Mass: {mass_low:.2e} Msun. High Mass: {mass_high:.2e} Msun.")
        if mass_low < 1e6: 
            mass_low = 1e6

        mask = (HaloMasses >= mass_low) & (HaloMasses <= mass_high)
        Nhalo = mask.sum()
        if Nhalo == 0:
            print("No halos found in the specified mass range.")
            return None

        nhalo = Nhalo / Volume

        print(f"Nhalo is {Nhalo} and nhalo is {nhalo}")
        mf = MassFunction(cosmo_model=cosmo, z=self.Redshift, Mmin=np.log10(mass_low)-1.0, Mmax=np.log10(mass_high)+1.0, dlog10m=0.01)

        mass = mf.m  # Msun/h (check units depending on settings)
        mask = (mass >= mass_low) & (mass <= mass_high)
        
        nhaloanalytic = np.trapz(mf.dndm[mask], mass[mask])

        norm = nhaloanalytic / nhalo
        print(f"n_halo (Sim) = {nhalo:.2e}, n_halo (Analytic) = {nhaloanalytic:.2e}, Normalization = {norm:.2f}")
        return norm

'''