import numpy as np
import os
from collections import defaultdict

from astropy.cosmology import FlatLambdaCDM
from GWFish.modules.detection import Network
from GWFish.modules.fishermatrix import compute_network_errors
from GWFish.modules.fishermatrix import compute_detector_fisher
from GWFish.modules.detection import Detector
import astropy.units as u

detector_lisa = Detector('LISA')
detector_et = Detector('ET')
detector_lgwa = Detector('LGWA')

network = Network(["LISA"])
cosmo = FlatLambdaCDM(H0=67.7, Om0=0.2592, Tcmb0=2.725)

from ..objects.merger import Merger, MergerCatalog
from ..core.constants import CS
from ..core.utils import find_snaps_around_merger

class MergerBuilder:

    def __init__(self, sim):
        self.sim = sim
    
    def build(self):

        sink_exclusion_list = self.sim.load_sink_exclusion_list()

        filename = os.path.join(self.sim.base, "sink_mergers.txt")
        seen = set()
        duplicates = 0

        mergers = []
        
        dtype = [
            ("Time", float),
            ("Redshift", float),

            ("ID1", int),
            ("ID2", int),

            ("Type1", int),
            ("Type2", int),

            ("M1", float),
            ("M2", float),

            ("P1x", float), ("P1y", float), ("P1z", float),
            ("P2x", float), ("P2y", float), ("P2z", float),

            ("V1x", float), ("V1y", float), ("V1z", float),
            ("V2x", float), ("V2y", float), ("V2z", float),
        ]

        data = np.genfromtxt(filename, names=True, dtype=dtype)

        for row in data:

            if row["ID1"] in sink_exclusion_list or row["ID2"] in sink_exclusion_list:
                continue

            key = (row["ID1"], row["ID2"])
            
            if key in seen:
                duplicates += 1
                continue

            seen.add(key)
            merger = Merger(time=row["Time"],  redshift=row["Redshift"], primary_id=row["ID1"], secondary_id= row["ID2"],  primary_type= row["Type1"], secondary_type= row["Type2"], primary_mass= row["M1"], secondary_mass= row["M2"], snapshot=find_snaps_around_merger(row["Redshift"], self.sim.snaps, self.sim.redshifts))
        
            mergers.append(merger)
        
        if duplicates:
            print(f"Removed {duplicates} duplicate merger records")
        
        print(f"We have a total of {len(mergers)} mergers.")

        self._calculate_snr(mergers)

        return MergerCatalog(mergers=mergers)
    
    @staticmethod
    def _calculate_snr(mergers):

        for merger in mergers:
            if merger.primary_type == 2 or merger.secondary_type == 2: 
                merger.snr_lisa = 0
                merger.snr_et = 0
                merger.snr_lgwa = 0
                continue

            params= {
                'mass_1': merger.primary_mass,
                'mass_2': merger.secondary_mass,
                'redshift': merger.redshift,
                'luminosity_distance': cosmo.luminosity_distance(merger.redshift).to(u.Mpc).value,
                #'theta_jn': 0.,
                'ra': 0.,
                'dec': -0.41,
                'phase': 0.,
                'psi': 1.6,
                'geocent_time': 1e9,
            }
            fisher, SNR_LISA = compute_detector_fisher(detector_lisa, params)
            fisher, SNR_ET = compute_detector_fisher(detector_et, params)
            fisger, SNR_LGWA = compute_detector_fisher(detector_lgwa, params)

            merger.snr_lisa = np.sqrt(SNR_LISA)
            merger.snr_et = np.sqrt(SNR_ET)
            merger.snr_lgwa = np.sqrt(SNR_LGWA)

            print(f"Calculated SNR for M1: {merger.primary_mass} and M2:{merger.secondary_mass}. SNR_LISA: {np.sqrt(SNR_LISA):.3f}, SNR_ET: {np.sqrt(SNR_ET):.3f}, SNR_LGWA: {np.sqrt(SNR_LGWA):.3f}")
