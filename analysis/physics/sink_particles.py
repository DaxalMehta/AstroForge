import os
import numpy as np
import pandas as pd
from astropy.cosmology import FlatLambdaCDM
from GWFish.modules.detection import Network
from GWFish.modules.fishermatrix import compute_network_errors
from GWFish.modules.fishermatrix import compute_detector_fisher
from GWFish.modules.detection import Detector
import astropy.units as u

detector = Detector('LISA')
network = Network(["LISA"])

from ..core import CS
from ..core import utils as UT

cosmo = FlatLambdaCDM(H0=67.7, Om0=0.2592, Tcmb0=2.725)

class SinkParticles:
    
    def __init__(self, sim):
        
        self.sim = sim
        
        self.base = sim.base
        self.feedback = sim.feedback
        self.array_length = sim.array_length

        
        self.reader = sim.reader

        self.full = None
        self.reduced = None
        self.merger_map = None
        self.popii = None
        self.popii_merger_map = None

        ## derived arrays
        self.ids = None
        self.M_init = None
        self.M_final = None
        self.M_diff = None
        self.Type = None
        self.FormationSnap = None
        self.FormationTime = None
        self.Status = None
    
    def load_full(self, force_create=True):
        
        if self.full is None:
            filepath = self.base + "sink_particle.pkl"
            if os.path.exists(filepath) and force_create is False:
                self.full = self.reader.read_pickle(os.path.join(self.base, "sink_particle.pkl"))["data"]
            else:
                self.full = self.sim.binary_reader.read_sink_info()
                #self.full = self.reader.read_pickle(os.path.join(self.base, "sink_particle.pkl"))["data"]
        
        return self.full
    
    def load_reduced(self):
        
        if self.reduced is None:
            filepath = self.base + "sink_particle_reduced.pkl"
            if os.path.exists(filepath):
                self.reduced = self.reader.read_pickle(os.path.join(self.base, "sink_particle_reduced.pkl"))
            else:
                self.build_reduced()
                self.reduced = self.reader.read_pickle(os.path.join(self.base, "sink_particle_reduced.pkl"))
                
            self._build_arrays()
            
        return self.reduced
    
    def load_popii(self):

        if self.popii is None:
            popii = self.reader.read_pickle(os.path.join(self.base, "popii_particles.pkl"))
            self.popii = popii["data"]
            print(popii["meta"])
        
        return self.popii

    def load_mergers(self, force_create=False):

        if self.merger_map is None:
            if os.path.exists(os.path.join(self.base, "merger_map.pkl")) and force_create is False:
                self.merger_map = self.reader.read_pickle(os.path.join(self.base, "merger_map.pkl"))
            else:
                builder = SinkBuilder(self.sim)
                self.merger_map = builder.create_merger_map()
        
        return self.merger_map

    def load_popii_mergers(self, force_create=False):

        if os.path.exists(os.path.join(self.base, "popii_merger_map.pkl")) and force_create is False:
            self.popii_merger_map = self.reader.read_pickle(os.path.join(self.base, "popii_merger_map.pkl"))
        else:
            builder = SinkBuilder(self.base, self.feedback, self.array_length)
            self.popii_merger_map = builder.create_popii_merger_map()
        
        return self.popii_merger_map

    def iter_sinks(self):

        self.load_full()

        for sid, sink in self.full.items():
            yield sid, sink
    
    def _build_arrays(self):
        
        data = self.reduced

        self.ids = np.array(list(data.keys()))
        
        meta = [data[i]["meta"] for i in self.ids]
        
        self.M_init = np.array([m["M_init"] for m in meta])
        self.M_final = np.array([m["M_final"] for m in meta])
        self.M_diff = self.M_final - self.M_init
        
        self.Type = np.array([m["Type"] for m in meta])
        self.FormationSnap = np.array([m["FormationSnap"] for m in meta])
        self.FormationTime = np.array([m["FormationTime"] for m in meta])
        self.Status = np.array([m["Status"] for m in meta])

    def build_reduced(self):
        builder = SinkBuilder(self.base, feedback=self.feedback, array_length=self.array_length)
        builder.create_sink_particles_reduced()
    
    def filter_by_type(self, t):
        
        if self.Type is None:
            self.load_reduced()
        
        return self.Type == t
         
    def filter_by_accretion(self, threshold):

        if self.M_diff is None:
            self.load_reduced()
        
        return self.M_diff > threshold
    
    def max_mass(self, t=0):
        
        mask = self.filter_by_type(t)
        
        return np.max(self.M_final[mask])
    
    def preliminary_statistics(self):

        self.load_reduced()

        print(f"There are {len(self.ids)} sink particles in the simulation.")
 
        for t in ["PopIII", "PopII", "MBH"]:
            
            mask = self.filter_by_type(t)
            print(f"There are {mask.sum()} {t} stars in the simulation.")
        

        mask = self.filter_by_accretion(self.M_init)

        print(f"There are {mask.sum()} BHs that doubled their initial mass")
        print("***********************************************************")
    
    def get_growth_tracks(self, types=None, mass_unit=True):

        data = self.load_full()

        tracks = {}

        for sid, sink in data.items():

            meta = sink["meta"]

            if types is not None and meta["Type"] not in types:
                continue

            evo = sink["evolution"]

            times = np.array(sorted(evo.keys()))

            stellar_mass = np.array([evo[t]["StellarMass"] for t in times])
            merger_mass = np.array([evo[t]["MergerMass"] for t in times])

            if mass_unit:
                stellar_mass *= 1e10 / CS.hubble_parameter
                merger_mass *= 1e10 / CS.hubble_parameter
            
            tracks[sid] = {
                "time": times,
                "redshift": 1 / times -1,
                "stellar_mass": stellar_mass,
                "merger_mass": merger_mass,
                "type": meta["Type"],
            }
        
        return tracks

    def get_merger_events(self):
        
        self.load_mergers()

        events = []
        
        for sid, m in self.merger_map.items():
            for t, m1, m2, mr, snr in zip(m["Time"], m["M1"], m["M2"], m["Mr"], m["SNR"]):
                events.append({
                    "sink_id": sid,
                    "time": t,
                    "redshift": 1/t - 1,
                    "M1": m1,
                    "M2": m2,
                    "Mr": mr,
                    "SNR": snr
                })
        return events

    def get_merger_arrays(self):

        events = self.get_merger_events()

        z  = np.array([e["redshift"] for e in events])
        M1 = np.array([e["M1"] for e in events])
        M2 = np.array([e["M2"] for e in events])
        Mr = np.array([e["Mr"] for e in events])
        snr = np.array([e["SNR"] for e in events])

        mask = M2 > M1
        M1[mask], M2[mask] = M2[mask], M1[mask]

        q = M2 / M1

        return z, M1, M2, Mr, q, snr

class SinkBuilder:
    def __init__(self, sim):
        
        self.base = sim.base
        self.reader = sim.reader
        self.writer = sim.writer
        self.binary_reader = sim.binary_reader
        filepath = self.base + "sink_particle.pkl"
        self.sink_particles = self.binary_reader.read_sink_info()
           
    def create_sink_particles_reduced(self):
        reduced_file = "sink_particle_reduced.pkl"

        valid_ids = set(self.sink_particles.keys())
        meta_cache = {sid: self.sink_particles[sid]["meta"] for sid in valid_ids}
        print(f"Number of sink particles in the simulation are {len(valid_ids)}")
        with open(f"{self.base}sink_id.txt", "w") as f:
            f.write("SinkID\n")
            f.write("\n".join(map(str, valid_ids)))
        
        start_snap, end_snap = UT.file_range(self.base, mode=1)
        
        sink_particle_reduced = {}
        
        for snap in range(start_snap, end_snap):
            self.binary_reader.read_sink_snap(snap)
            time = self.binary_reader.time
            SinkID = self.binary_reader.extract_data("ID")
            SinkMass = self.binary_reader.extract_data("StellarMass")
            SinkType = self.binary_reader.extract_data("Type")
            SinkPos = self.binary_reader.extract_data("Pos")
            SinkVel = self.binary_reader.extract_data("Vel")

            for i, sid in enumerate(SinkID):

                if sid not in valid_ids:
                    continue
                mass, pos, vel, typ = SinkMass[i], SinkPos[i], SinkVel[i], SinkType[i]
                meta = meta_cache[sid]

                sink = sink_particle_reduced.setdefault(sid, {
                    "meta": {
                        "FormationSnap": snap,
                        "M_init": mass,
                        "M_final": mass,
                        "StellarLifeTime": meta["StellarLifeTime"],
                        "FormationTime": meta["FormationTime"],
                        "Type": meta["Type"],
                        "SNeType": meta["SNeType"],
                        "Status": meta["Status"],
                    },
                    "evolution": {}
                })

                sink["meta"]["M_final"] = mass
                
                sink["evolution"][time] = {
                    "Time": time,
                    "Pos": pos,
                    "Vel": vel,
                    "StellarMass": mass,
                    "Type": typ,
                }
                
        for sid in valid_ids:
            if sid not in sink_particle_reduced:

                sink = self.sink_particles[sid]
                meta = sink["meta"]
                evolution = sink["evolution"]
                times = np.array(sorted(evolution.keys()))
                init_time = times[0]
                final_time = times[-1]
                sink_particle_reduced[sid] = {
                    "meta": {
                        "FormationSnap": None,
                        "M_init": evolution[init_time]["StellarMass"],
                        "M_final": evolution[final_time]["StellarMass"],
                        "StellarLifeTime": meta["StellarLifeTime"],
                        "FormationTime": meta["FormationTime"],
                        "Type": meta["Type"], # Can be PopII, PopIII, or MBH
                        "SNeType": meta["SNeType"], # Can be No SNe, TypeII SNe, DCBH, PISN
                        "Status": meta["Status"],
                    },
                    "evolution": {}
                }
            if sink_particle_reduced[sid]["meta"]["M_final"] == None:
                evolution = self.sink_particles[sid]["evolution"]
                final_time = max(evolution)
                sink_particle_reduced[sid]["meta"]["M_final"] = evolution[final_time]["StellarMass"]
         
        self.writer.write_pickle(reduced_file, sink_particle_reduced)

    def create_merger_map(self):
        
        file = os.path.join(self.base, "merger_map.pkl")

        merger_map = {}
        
        filename = os.path.join(self.base, "preliminary_info.txt")
        snaps, times, redshifts = np.genfromtxt(filename, usecols=(0, 1, 2), skip_header=1, unpack=True)
        
        for ID, sink in self.sink_particles.items():
            
            evo = sink["evolution"]
            
            times = np.array(sorted(evo.keys()))
            
            merger_mass = np.array([evo[t]["MergerMass"] for t in times])
            stellar_mass = np.array([evo[t]["StellarMass"] for t in times])
            
            #merger increments
            deltas = np.diff(merger_mass)
            mask = deltas > 0
            if not np.any(mask):
                continue

            t_merge = times[:-1][mask]
            z_merge = 1 / t_merge - 1
            z_snaps = []
            for z in z_merge:
                z_snaps.append(int(UT.find_snaps_around_merger(z, snaps, redshifts)[0]))
            
            z_snaps = np.array(z_snaps)
            m2 = deltas[mask]
            m1 = stellar_mass[:-1][mask]

            snr = []

            for primary, secondary, z in zip(m1, m2, z_merge):
                params= {
                    'mass_1': primary * 1e10 / CS.hubble_parameter,
                    'mass_2': secondary * 1e10 / CS.hubble_parameter,
                    'redshift': z,
                    'luminosity_distance': cosmo.luminosity_distance(z).to(u.Mpc).value,
                    #'theta_jn': 0.,
                    'ra': 0.,
                    'dec': -0.41,
                    'phase': 0.,
                    'psi': 1.6,
                    'geocent_time': 1e9,
                }
                fisher, SNR = compute_detector_fisher(detector, params)
                #parameters = pd.DataFrame.from_dict({k:v*np.array([1.]) for k, v in param_dict.items()})
                #detected, SNR, errors, sky_localization = compute_network_errors(network, parameters, waveform_model='TaylorF2')
                snr.append(np.sqrt(SNR))
                print(SNR)
            
            snr = np.array(snr)
            
            merger_map[ID] = {
                "M1": m1,
                "M2": m2,
                "Mr": m1 + m2,
                "Time": t_merge,
                "Snap": z_snaps,
                "Mi": stellar_mass[0],
                "Mf": stellar_mass[-1],
                "SNR": snr
            } 
        
        self.writer.write_pickle(file, merger_map)

        return merger_map

    def create_popii_merger_map(self):

        file = os.path.join(self.base, "popii_merger_map.pkl")

        merger_map = {}
        popii_particles = self.reader.read_pickle(os.path.join(self.base, "popii_particles.pkl"))["data"]

        for ID, popii in popii_particles.items():
            evo = popii["evolution"]
            
            times = np.array(sorted(evo.keys()))
            
            merger_mass = np.array([evo[t]["MergerMass"] for t in times])
            stellar_mass = np.array([evo[t]["StellarMass"] for t in times])
            
            #merger increments
            deltas = np.diff(merger_mass)
            mask = deltas > 0
            if not np.any(mask):
                continue

            t_merge = times[:-1][mask]
            z_merge = 1 / t_merge - 1

            m2 = deltas[mask]
            m1 = stellar_mass[:-1][mask]
            
            merger_map[ID] = {
                "M1": m1,
                "M2": m2,
                "Mr": m1 + m2,
                "Time": t_merge,
                "Mi": stellar_mass[0],
                "Mf": stellar_mass[-1]
            } 
        
        self.Writer.write_pickle(file, merger_map)

        return merger_map

class SinkAnalysis:
    def __init__(self, sim):
        
        self.base = sim.base
        self.reader = sim.reader
        self.binary_reader = sim.binary_reader

        self.sink_particles = self.reader.read_pickle("sink_particle.pkl")
    
    def bh_mass_function(self, box_volume, IDs=None):
        if IDs == None:
            IDs = np.array(list(self.sink_particles.keys()))
        MBHMass = []
        MBHAccretion = []
        MBHLuminosity = []
        for i, id in enumerate(IDs):
            if id == -1: continue
            if id not in self.sink_particles: continue
            evolution = self.sink_particles[id]["evolution"]
            times = np.array(sorted(evolution.keys()))
            Redshift = 1 / times - 1
            EvolutionTimes = Time = cosmo.age(Redshift).to("Myr").value
            mass = np.array([evolution[t]["StellarMass"] for t in times]) * 1e10 / CS.hubble_parameter
            if len(mass) < 2: continue
            if mass[-1] < 1e5: continue
            MBHMass.append(mass[-1])
            if mass[-1] > mass[-2]: MBHAccretion.append(mass[-1])
            massdiff = (mass[-1] - mass[-2]) * CS.Msun_to_g
            timediff = (EvolutionTimes[-1] - EvolutionTimes[-2]) * 1e6 * CS.yr_to_s
            AccretionRate = massdiff/ timediff

            Luminosity = 0.1 * AccretionRate * pow(CS.c_in_cgs, 2)
            if Luminosity > 1e42:
                MBHLuminosity.append(mass[-1])
        
        log_min = 5.25 #np.min(np.log10(MBHMass))
        log_max = np.max(np.log10(MBHMass))
        bin_width = 0.50
        log_centers = np.arange(log_min, log_max + bin_width, bin_width)
        log_edges = log_centers - bin_width / 2
        log_edges = np.append(log_edges, log_edges[-1] + bin_width)
        bin_edges = 10**log_edges
        bin_centers = 10**log_centers
        counts1, _ = np.histogram(MBHMass, bins=bin_edges)
        counts2, _ = np.histogram(MBHAccretion, bins=bin_edges)
        counts3, _ = np.histogram(MBHLuminosity, bins=bin_edges)

        bin_wdith_dex = bin_width

        dN_dM1 = counts1 / (box_volume * bin_wdith_dex) 
        dN_dM2 = counts2 / (box_volume * bin_wdith_dex)
        dN_dM3 = counts3 / (box_volume * bin_wdith_dex)

        return dN_dM1, dN_dM2, dN_dM3, bin_centers

    def luminosity_function(self):
        sink_id = np.array(list(self.sink_particles.keys()))
        MBHMass = []
        Lbol = []
        for i, id in enumerate(sink_id):
            if self.sink_particles[id]["meta"]["Type"] != "MBH": continue
            evolution = self.sink_particles[id]["evolution"]
            times = np.array(sorted(evolution.keys()))
            Redshift = 1 / times - 1
            EvolutionTimes = Time = cosmo.age(Redshift).to("Myr").value
            mass = np.array([evolution[t]["StellarMass"] for t in times]) * 1e10 / CS.hubble_parameter
            if len(mass) < 2: continue
            if mass[-1] < 1e5: continue

            massdiff = (mass[-1] - mass[-2]) * CS.Msun_to_g
            timediff = (EvolutionTimes[-1] - EvolutionTimes[-2]) * 1e6 * CS.yr_to_s
            AccretionRate = massdiff / timediff

            Luminosity = 0.1 * AccretionRate * pow(CS.c_in_cgs, 2)
            if Luminosity > 1e42:
                MBHMass.append(mass[-1])
                Lbol.append(Luminosity)
        
        MBHMass, Lbol = np.array(MBHMass), np.array(Lbol)

        return MBHMass, Lbol