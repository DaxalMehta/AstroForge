import numpy as np
import os
from dataclasses import dataclass

from ..core.enums import SinkType
from ..core.utils import unwrap
from ..objects.metadata import ParticleDatabase
from ..objects.particle import ParticleHistory, ParticleHistoryDatabase
from ..inout.binary_reader import BinaryReader, BinaryReaderSnapshot

@dataclass
class ParticleBuildState:
    files_processed: list[str]
    last_dump_time: float
    start_times: list[float]
    end_times: list[float]

    def to_dict(self):
        """
        Convert the particle build state class to a dictionary for saving.
        """

        return {
            "files_processed": self.files_processed,
            "last_dump_time": self.last_dump_time,
            "start_times": self.start_times,
            "end_times": self.end_times
        }
    
    @classmethod
    def from_dict(cls, data):
        """
        Construct a ParticleBuildState from a dictionary produced 
        by 'to_dict()'
        """
        return cls(**data)

class ParticleBuilder:

    def __init__(self, sim):
        self.sim = sim
        self.base = sim.base
        self.binary_reader = BinaryReader(sim)

        self._seen = set()
        self._ids = []
        self._formation_times = []
        self._initial_masses = []
        self._particle_types = []
        self._stellar_lifetimes = []

        self._history = ParticleHistoryDatabase()
    
    def _find_files(self):

        sink_base = f"{self.base}sink_particle_info/"

        return sorted(os.path.join(sink_base, file) for file in os.listdir(sink_base))

    def _load_build_state(self):
        filename = "particle_build_state.pkl"

        if os.path.exists(os.path.join(self.base, filename)):
            data =  self.sim.reader.read_pickle(filename)
            return ParticleBuildState.from_dict(data)
             
        return None
    
    def _save_build_state(self, files_processed, last_dump_time, start_times, end_times):
        state = ParticleBuildState(
            files_processed=files_processed,
            last_dump_time=last_dump_time,
            start_times=start_times,
            end_times=end_times
        )
        self.sim.writer.write_pickle("particle_build_state.pkl", state.to_dict())

    def _update_database(self, record):

        if record.id in self._seen:
            return

        self._seen.add(record.id)

        self._ids.append(record.id)
        self._initial_masses.append(record.mass)
        self._particle_types.append(record.particle_type)

        self._formation_times.append(record.formation_time)
        self._stellar_lifetimes.append(record.stellar_lifetime)

    def _update_history(self, record):
        if record.particle_type == 3:
            if record.id not in self._history.black_holes:
                self._history.black_holes[record.id] = ParticleHistory(record.id, record.particle_type)
            self._history.black_holes[record.id].append(record)
        else:
            if record.id not in self._history.stellar_particles:
                self._history.stellar_particles[record.id] = ParticleHistory(record.id, record.particle_type)
            self._history.stellar_particles[record.id].append(record)

    def _rollback_to_time(self, restart_time):

        """
        Removes all particle information that was created after a simulation restart. This allows the new binary files to rebuild the correct history.
        """
        # ----------------------
        # Particle histories
        # ----------------------

        def truncate(history_dict):
            remove = []
            for pid, history in history_dict.items():
                mask = np.asarray(history.dump_times) < restart_time

                if np.any(mask):
                    history.dump_times = list(np.asarray(history.dump_times)[mask])
                    history.masses = list(np.asarray(history.masses)[mask])
                    history.merger_masses = list(np.asarray(history.merger_masses)[mask])
                    history.positions = list(np.asarray(history.positions)[mask])
                    history.velocities = list(np.asarray(history.velocities)[mask])
                
                else:
                    remove.append(pid)
            
            for pid in remove:
                del history_dict[pid]

        truncate(self._history.black_holes)
        truncate(self._history.stellar_particles)

        # ------------------------
        # Metadata
        # ------------------------

        keep = np.asarray(self._formation_times) <= restart_time

        self._ids = list(np.asarray(self._ids)[keep])
        self._formation_times = list(np.asarray(self._formation_times)[keep])
        self._initial_masses = list(np.asarray(self._initial_masses)[keep])
        self._particle_types = list(np.asarray(self._particle_types)[keep])
        self._stellar_lifetimes = list(np.asarray(self._stellar_lifetimes)[keep])

        self._seen = set(self._ids)

    def build(self):

        files = self._find_files()

        state = self._load_build_state()

        if state is None:
            files_to_process = files
            last_dump_time = -np.inf
            start_times = []
            end_times = []

        else:
            new_files = [
                file for file in files if file not in state.files_processed
            ]

            if len(new_files) == 0:
                print("Found no new files to process. Loading existing pickle files.")

                data = self.sim.reader.read_pickle("metadata.pkl")
                metadata = ParticleDatabase.from_dict(data)

                history = ParticleHistoryDatabase()
                black_hole_data = self.sim.reader.read_pickle("black_holes.pkl")
                history.black_holes = {
                    int(pid): ParticleHistory.from_dict(history_data)
                    for pid, history_data in black_hole_data.items()
                }
                stellar_particles_data = self.sim.reader.read_pickle("stellar_particles.pkl")
                history.stellar_particles = {
                    int(pid): ParticleHistory.from_dict(history_data)
                    for pid, history_data in stellar_particles_data.items()
                }

                return metadata, history
            
            else:
                print("Found new files to process. Running incremental updates.")
                
                data = self.sim.reader.read_pickle("metadata.pkl")
                metadata = ParticleDatabase.from_dict(data)

                self._ids = list(metadata.ids)
                self._initial_masses = list(metadata.initial_masses)
                self._particle_types = list(metadata.particle_types)
                self._formation_times = list(metadata.formation_times)
                self._stellar_lifetimes = list(metadata.stellar_lifetimes)

                self._seen = set(metadata.ids)

                history = ParticleHistoryDatabase()
                black_hole_data = self.sim.reader.read_pickle("black_holes.pkl")
                history.black_holes = {
                    int(pid): ParticleHistory.from_dict(history_data)
                    for pid, history_data in black_hole_data.items()
                }
                stellar_particles_data = self.sim.reader.read_pickle("stellar_particles.pkl")
                history.stellar_particles = {
                    int(pid): ParticleHistory.from_dict(history_data)
                    for pid, history_data in stellar_particles_data.items()
                }

                files_to_process = new_files

                last_dump_time = state.last_dump_time
                start_times = state.start_times
                end_times = state.end_times


        for file in files_to_process:
            sink_dumps = self.binary_reader.read_sink_info(file)

            if len(sink_dumps) == 0:
                continue

            start_times.append(sink_dumps[0].time)
            end_times.append(sink_dumps[-1].time)

            first_time = sink_dumps[0].time
            if first_time < last_dump_time:
                print(f"Restart detected at t = {first_time:.6f}")
                self._rollback_to_time(first_time)

                last_dump_time = first_time

            for dump in sink_dumps:

                if dump.time < last_dump_time:
                    continue

                last_dump_time = dump.time

                for record in dump.records:
                    self._update_database(record)
                    self._update_history(record)
                   
        metadata = ParticleDatabase(
            ids=self._ids,
            initial_masses=self._initial_masses,
            particle_types=self._particle_types,
            formation_times=self._formation_times,
            stellar_lifetimes=self._stellar_lifetimes
        )
        
        self.sim.writer.write_pickle("metadata.pkl", metadata.to_dict())
        self.sim.writer.write_pickle("black_holes.pkl", 
                                     {
                                        particle_id: history.to_dict()
                                        for particle_id, history in self._history.black_holes.items()
                                     })
        self.sim.writer.write_pickle("stellar_particles.pkl", 
                                     {
                                        particle_id: history.to_dict()
                                        for particle_id, history in self._history.stellar_particles.items()
                                     })

        self._save_build_state(files_processed=files, last_dump_time=last_dump_time, start_times=start_times, end_times=end_times)
                    
        return metadata, self._history

    def complete_particle_database(self):
        """
        Give Final Masses to the particle database
        """

        end_snap = self.sim.end_snap - 1

        data = self.sim.reader.read_pickle("metadata.pkl")
        metadata = ParticleDatabase.from_dict(data)
        
        sink_snap = BinaryReaderSnapshot(self.sim).read_sink_snap(end_snap)

        mass_lookup = dict(zip(sink_snap.ids, sink_snap.masses))

        final_masses = np.array(
            [mass_lookup.get(pid, np.nan) for pid in metadata.ids]
        )

        metadata.final_masses = final_masses
        
        self.sim.writer.write_pickle("metadata.pkl", metadata.to_dict())
        




