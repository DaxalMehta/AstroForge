import numpy as np
import os
import pickle
import struct
from ..core.utils import Type, SNeType
from ..core.enums import SinkType
from .binary_formats import sink_info
from .binary_records import SinkInfoRecord, SinkDump, SinkSnapshot, SinkSnapshotRecord

class BinaryReader:
    
    def __init__(self, sim):
        self.sim = sim
        self.base = sim.base
        self.feedback = sim.feedback
        
        self.format = sink_info(feedback=sim.feedback, array_length=sim.array_length)

    def _read_record(self, file):

        raw = file.read(self.format.size)

        if len(raw) != self.format.size:
            raise EOFError
        
        return struct.unpack(self.format.format, raw)
    
    def _create_record(self, dump_time, data):

        return SinkInfoRecord(
            dump_time=dump_time, 

            position = np.asarray(data[0:3]),
            velocity=np.asarray(data[3:6]),
            acceleration=np.asarray(data[6:9]),

            mass=data[9],

            formation_time=data[12],

            id=data[13],

            particle_type=data[17],

            stellar_lifetime=data[19],

            merger_mass=data[24] if self.feedback else data[22],
        )

    def read_sink_info(self, filename):
        print(f"Opening file {filename}")
        sink_dumps = []
        
        with open(filename, "rb") as f:
            while True:
                chunk = f.read(8)
                if not chunk or len(chunk) < 8:
                    break
                time = struct.unpack("d", chunk)[0]
                if (time < 0) or (time > 1):
                    break
                num_sinks = struct.unpack("i", f.read(4))[0]
                if (num_sinks < 0):# or (num_sinks > 50000):
                    break
                print(f"Time of the dump: {time:.4f}. Number of sink particles: {num_sinks}")
                records = []
                for _ in range(num_sinks):
                    try:
                        data = self._read_record(f)
                    except EOFError:
                        print(f"Incomplete sink entry at time {time}, aborting for this timestep.")
                        break

                    record = self._create_record(time, data)

                    records.append(record)
        

                sink_dumps.append(SinkDump(
                    time = time,
                    records=records,
                ))
        
        
        return sink_dumps

class BinaryReaderSnapshot:

    def __init__(self, sim):
        self.sim = sim
        self.base = sim.base
        self.feedback = sim.feedback
        
        self.format = sink_info(feedback=sim.feedback, array_length=sim.array_length)
    
    def _read_record(self, file):

        raw = file.read(self.format.size)

        if len(raw) != self.format.size:
            raise EOFError
        
        return struct.unpack(self.format.format, raw)
    
    def _create_record(self, data):

        return SinkSnapshotRecord(

            position = np.asarray(data[0:3]),
            velocity=np.asarray(data[3:6]),
            acceleration=np.asarray(data[6:9]),

            mass=data[9],

            formation_time=data[12],

            id=data[13],

            particle_type=data[17],

            stellar_lifetime=data[19],

            merger_mass=data[24] if self.feedback else data[22],

            number_of_sne=data[26] if self.feedback else data[24]  
        )

    def read_sink_snap(self, snap):
        filename = os.path.join(self.base, f"sink_snap_{snap:03d}")
        records = []
        with open(filename, "rb") as f:
            time = struct.unpack("d", f.read(8))[0]
            redshift = 1 / time - 1
            num_sinks = struct.unpack("i", f.read(4))[0]

            for _ in range(num_sinks):
                data = self._read_record(f)

                record = self._create_record(data)

                records.append(record)
        
        return SinkSnapshot(time=time, redshift=redshift, sink_particles=records)


'''
    def read_sink_info(self):
        
        sink_base = f"{self.base}sink_particle_info/"
        files = sorted(os.listdir(sink_base))
        
        pickle_file = f"{self.base}sink_particle.pkl"
        pickle_file_popii = f"{self.base}popii_particles.pkl"

        if os.path.exists(pickle_file):
            with open(pickle_file, "rb") as f:
                saved = pickle.load(f)
            
            saved_files = saved["meta"]["files_used"]

            new_files = [f for f in files if f not in saved_files]

            if len(new_files) == 0:
                print("No new files. Loading existing pickle.")

                sink_particles = saved["data"]
                #with open(pickle_file, "rb") as f:
                #    sink_particles = pickle.load(f)["data"]
                #with open(pickle_file_popii, "rb") as f:
                #    popii_particles = pickle.load(f)["data"]
                
                return sink_particles#, popii_particles
            
            else:
                print(f"New files detected: {len(new_files)}. Running incremental update.")

                #load existing data
                sink_particles = saved["data"]
                with open(pickle_file_popii, "rb") as f:
                    popii_particles = pickle.load(f)["data"]
                
                files_to_process = new_files
        
        else: # No pickle exists
            print("No pickle found. Running full build.")

            sink_particles = {}
            popii_particles = {}
            files_to_process = files

        if "saved" in locals():
            last_time = saved["meta"].get("last_time", -np.inf)
        else:
            last_time = -np.inf

        for fname in files_to_process:
            filename = f"{sink_base}{fname}"
            print(f"Opening file {filename}")

            with open(filename, "rb") as f:
                while True:
                    chunk = f.read(8)
                    if not chunk or len(chunk) < 8:
                        break
                    time = struct.unpack("d", chunk)[0]
                    if (time < 0) or (time > 1):
                        break
                    num_sinks = struct.unpack("i", f.read(4))[0]
                    if (num_sinks < 0):# or (num_sinks > 50000):
                        break
                    if time <= last_time:
                        f.seek(self.struct_size * num_sinks, 1)
                        continue
                    print(f"Time of the dump: {time}. Number of sink particles: {num_sinks}")

                    for _ in range(num_sinks):
                        try:
                            raw = f.read(self.struct_size)
                            if len(raw) < self.struct_size:
                                raise struct.error("Incomplete data.")
                                break
                            data = struct.unpack(self.struct_format, raw)
                        except struct.error:
                            print(f"Incomplete sink entry at time {time}, aborting for this timestep.")
                            break
                        
                        sink_id = data[13]
                        entry = {
                            "Time": time,
                            "Pos": [data[0:3]],
                            "Vel": [data[3:6]],
                            "Type": data[17],
                            "StellarMass": data[23] if self.feedback else data[21],
                            "MergerMass": data[24] if self.feedback else data[22],
                        }                        
                        if entry["Type"] in [0, 3]: 
                            if sink_id not in sink_particles:
                                print(f"Adding sink_particle ID {sink_id}")
                                if sink_id < 0: break
                                sink_particles[sink_id] = {
                                    "meta": {
                                        "StellarLifeTime": data[19],
                                        "FormationTime": data[12],
                                        "Type": Type(data[17]), # Can be PopII, PopIII, or MBH
                                        "SNeType": SNeType(data[9]) if data[17] == 0 else Type(data[17]), # Can be No SNe, TypeII SNe, DCBH, PISN
                                        "Status": None,
                                        
                                    },
                                    "evolution": {}
                                }
                            sink_particles[sink_id]["evolution"][time] = entry
                        if entry["Type"] == 2:
                            if sink_id not in popii_particles:
                                print(f"Adding PopII particle ID {sink_id}")
                                if sink_id < 0: break
                                popii_particles[sink_id] = {
                                    "meta": {
                                        "StellarLifeTime": data[19],
                                        "FormationTime": data[12],
                                        "Type": Type(data[17]),
                                        "Status": None,
                                    },
                                    "evolution": {}
                                }
                            popii_particles[sink_id]["evolution"][time] = entry

                    last_time = time             
        all_times = []
        for sink_id in sink_particles:
            sink_particles[sink_id]["evolution"] = dict(sorted(sink_particles[sink_id]["evolution"].items()))
            all_times.extend(sink_particles[sink_id]["evolution"].keys())

        for sink_id in popii_particles:
            popii_particles[sink_id]["evolution"] = dict(sorted(popii_particles[sink_id]["evolution"].items()))
        
        final_time = max(all_times)
        for sink_id in sink_particles:
            meta = sink_particles[sink_id]["meta"]
            evolution = sink_particles[sink_id]["evolution"]
            last_present = final_time in evolution
                
            if meta["Type"] == "PopIII":
                if meta["SNeType"] == "PISN":
                    meta["Status"] = "In Simulation" if last_present else "Deleted"
                else:
                    meta["Status"] = "In Simulation" if last_present else "Merged"
            elif meta["Type"] == "MBH":
                meta["Status"] = "In Simulation" if last_present else "Merged"
        
        for sink_id in popii_particles:
            meta = popii_particles[sink_id]["meta"]
            evolution = popii_particles[sink_id]["evolution"]
            last_present = final_time in evolution

            meta["Status"] = "In Simulation" if last_present else "Merged"

        pickle_file = f"{self.base}sink_particle.pkl"
        pickle_file_popii = f"{self.base}popii_particles.pkl"

        
        metadata = {
            "files_used": files,
            "last_time": final_time
        }


        with open(pickle_file, "wb") as f:
            #pickle.dump(sink_particles, f)
            pickle.dump({"data": sink_particles, "meta": metadata}, f)
        
        with open(pickle_file_popii, "wb") as f:
            #pickle.dump(popii_particles, f)
            pickle.dump({"data": popii_particles, "meta": metadata}, f)
        
        print("sink_particle.pkl and popii_partickles.pkl file created.")

        return sink_particles

    def read_sink_snap(self, snap, log=True):
        filename = f"{self.base}sink_snap_{snap:03d}"
        with open(filename, "rb") as f:
            time = struct.unpack("d", f.read(8))[0]
            Redshift = 1 / time - 1
            num_sinks = struct.unpack("i", f.read(4))[0]
            if log:
                print(f"Time of snapshot: {time}. Redshift of the snapshot: {Redshift}. Number of sink particles: {num_sinks}")

            sink_particles = []
            for _ in range(num_sinks):
                data = struct.unpack(self.struct_format, f.read(self.struct_size))

                sink_particles.append({
                    "Pos": data[0:3],
                    "Vel": data[3:6],
                    "Mass": data[9],
                    "ID": data[13],
                    "FormationTime": data[12],
                    "Type": data[17],
                    "StellarMass" : data[23] if self.feedback else data[21],
                    "MergerMass": data[24] if self.feedback else data[22],
                    "N_SNe": data[26] if self.feedback else data[24]
                })
        self.time = time
        self.num_sinks = num_sinks
        self.sink_particles = sink_particles
    
    def extract_data(self, field):
        return np.array([sink[field] for sink in self.sink_particles])   


class BinaryReader_IsolatedGalaxy:
    
    def __init__(self, sim):
        self.sim = sim
        self.base = sim.base
        self.feedback = sim.feedback
    
        if self.feedback:
            main = (
                "3d 3d 3d d d d d d d d i i i i d d d d i d d d d"
            )
            self.struct_format = main 
        else:
            main = (
                "3d 3d 3d d d d d d d d i i i i d d d d"
            )
            self.struct_format = main
            
        self.struct_size = struct.calcsize(self.struct_format)

    def read_sink_info(self):
        
        sink_base = f"{self.base}sink_particle_info/"
        files = sorted(os.listdir(sink_base))
        
        pickle_file = f"{self.base}sink_particle.pkl"
        pickle_file_popii = f"{self.base}popii_particles.pkl"

        if os.path.exists(pickle_file):
            with open(pickle_file, "rb") as f:
                saved = pickle.load(f)
            
            saved_files = saved["meta"]["files_used"]

            new_files = [f for f in files if f not in saved_files]

            if len(new_files) == 0:
                print("No new files. Loading existing pickle.")

                sink_particles = saved["data"]
                #with open(pickle_file, "rb") as f:
                #    sink_particles = pickle.load(f)["data"]
                #with open(pickle_file_popii, "rb") as f:
                #    popii_particles = pickle.load(f)["data"]
                
                return sink_particles#, popii_particles
            
            else:
                print(f"New files detected: {len(new_files)}. Running incremental update.")

                #load existing data
                sink_particles = saved["data"]
                with open(pickle_file_popii, "rb") as f:
                    popii_particles = pickle.load(f)["data"]
                
                files_to_process = new_files
        
        else: # No pickle exists
            print("No pickle found. Running full build.")

            sink_particles = {}
            popii_particles = {}
            files_to_process = files

        if "saved" in locals():
            last_time = saved["meta"].get("last_time", -np.inf)
        else:
            last_time = -np.inf

        for fname in files_to_process:
            filename = f"{sink_base}{fname}"
            print(f"Opening file {filename}")

            with open(filename, "rb") as f:
                while True:
                    chunk = f.read(8)
                    if not chunk or len(chunk) < 8:
                        break
                    time = struct.unpack("d", chunk)[0]
                    if (time < 0) or (time > 1):
                        break
                    num_sinks = struct.unpack("i", f.read(4))[0]
                    if (num_sinks < 0):# or (num_sinks > 50000):
                        break
                    if time <= last_time:
                        f.seek(self.struct_size * num_sinks, 1)
                        continue
                    print(f"Time of the dump: {time}. Number of sink particles: {num_sinks}")

                    for _ in range(num_sinks):
                        try:
                            raw = f.read(self.struct_size)
                            if len(raw) < self.struct_size:
                                raise struct.error("Incomplete data.")
                                break
                            data = struct.unpack(self.struct_format, raw)
                        except struct.error:
                            print(f"Incomplete sink entry at time {time}, aborting for this timestep.")
                            break
                        
                        sink_id = data[13]
                        entry = {
                            "Time": time,
                            "Pos": [data[0:3]],
                            "Vel": [data[3:6]],
                            "Type": data[17],
                            "StellarMass": data[23] if self.feedback else data[21],
                            "MergerMass": data[24] if self.feedback else data[22],
                        }                        
                        if entry["Type"] in [0, 3]: 
                            if sink_id not in sink_particles:
                                print(f"Adding sink_particle ID {sink_id}")
                                if sink_id < 0: break
                                sink_particles[sink_id] = {
                                    "meta": {
                                        "StellarLifeTime": data[19],
                                        "FormationTime": data[12],
                                        "Type": Type(data[17]), # Can be PopII, PopIII, or MBH
                                        "SNeType": SNeType(data[9]) if data[17] == 0 else Type(data[17]), # Can be No SNe, TypeII SNe, DCBH, PISN
                                        "Status": None,
                                        
                                    },
                                    "evolution": {}
                                }
                            sink_particles[sink_id]["evolution"][time] = entry
                        if entry["Type"] == 2:
                            if sink_id not in popii_particles:
                                print(f"Adding PopII particle ID {sink_id}")
                                if sink_id < 0: break
                                popii_particles[sink_id] = {
                                    "meta": {
                                        "StellarLifeTime": data[19],
                                        "FormationTime": data[12],
                                        "Type": Type(data[17]),
                                        "Status": None,
                                    },
                                    "evolution": {}
                                }
                            popii_particles[sink_id]["evolution"][time] = entry

                    last_time = time             
        all_times = []
        for sink_id in sink_particles:
            sink_particles[sink_id]["evolution"] = dict(sorted(sink_particles[sink_id]["evolution"].items()))
            all_times.extend(sink_particles[sink_id]["evolution"].keys())

        for sink_id in popii_particles:
            popii_particles[sink_id]["evolution"] = dict(sorted(popii_particles[sink_id]["evolution"].items()))
        
        final_time = max(all_times)
        for sink_id in sink_particles:
            meta = sink_particles[sink_id]["meta"]
            evolution = sink_particles[sink_id]["evolution"]
            last_present = final_time in evolution
                
            if meta["Type"] == "PopIII":
                if meta["SNeType"] == "PISN":
                    meta["Status"] = "In Simulation" if last_present else "Deleted"
                else:
                    meta["Status"] = "In Simulation" if last_present else "Merged"
            elif meta["Type"] == "MBH":
                meta["Status"] = "In Simulation" if last_present else "Merged"
        
        for sink_id in popii_particles:
            meta = popii_particles[sink_id]["meta"]
            evolution = popii_particles[sink_id]["evolution"]
            last_present = final_time in evolution

            meta["Status"] = "In Simulation" if last_present else "Merged"

        pickle_file = f"{self.base}sink_particle.pkl"
        pickle_file_popii = f"{self.base}popii_particles.pkl"

        
        metadata = {
            "files_used": files,
            "last_time": final_time
        }


        with open(pickle_file, "wb") as f:
            #pickle.dump(sink_particles, f)
            pickle.dump({"data": sink_particles, "meta": metadata}, f)
        
        with open(pickle_file_popii, "wb") as f:
            #pickle.dump(popii_particles, f)
            pickle.dump({"data": popii_particles, "meta": metadata}, f)
        
        print("sink_particle.pkl and popii_partickles.pkl file created.")

        return sink_particles

    def read_sink_snap(self, snap, log=True):
        filename = f"{self.base}sink_snap_{snap:03d}"
        with open(filename, "rb") as f:
            time = struct.unpack("d", f.read(8))[0]
            Redshift = 1 / time - 1
            num_sinks = struct.unpack("i", f.read(4))[0]
            if log:
                print(f"Time of snapshot: {time}. Redshift of the snapshot: {Redshift}. Number of sink particles: {num_sinks}")

            sink_particles = []
            for _ in range(num_sinks):
                data = struct.unpack(self.struct_format, f.read(self.struct_size))

                sink_particles.append({
                    "Pos": data[0:3],
                    "Vel": data[3:6],
                    "Mass": data[9],
                    "ID": data[14],
                })
        self.time = time
        self.num_sinks = num_sinks
        self.sink_particles = sink_particles
    
    def extract_data(self, field):
        return np.array([sink[field] for sink in self.sink_particles])   
''' 