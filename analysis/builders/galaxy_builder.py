import os
import numpy as np

from ..core.constants import CS
from ..objects.galaxy import Galaxy, GalaxyCatalog

class CatalogBuilder:
    
    def __init__(self, sim):
        
        self.sim = sim
        self.base = sim.base

    def filename(self, snap):
        return f"Galaxies_{snap:03d}.pkl"
    
    def build_catalog(self, snap):
        outfile = self.filename(snap)

        if os.path.exists(os.path.join(self.base, outfile)):
            print(f"Galaxy catalog already exists fro snap {snap}")
            return
        
        snapshot = self.sim.snapshot(snap)
        
        groups = snapshot.groups
        subhalos = groups.subhalos
        
        if groups.header.Ngroups_ThisFile == 0:
            print(f"There are no halos in the given snapshot. Skipping it.")
            return


        fof_indexes = subhalos.group_index

        position = subhalos.position
        velocity = subhalos.velocity
        center_of_mass = subhalos.center_of_mass

        catalog_mass = subhalos.mass * 1e10 / CS.hubble_parameter
        catalog_mass_type = subhalos.mass_type * 1e10 / CS.hubble_parameter

        halfmass_radius = subhalos.halfmass_radius
        halfmass_radius_type = subhalos.halfmass_radius_type
        aperture_radius = 2 * halfmass_radius_type[:, 4]

        vmax = subhalos.vmax
        vmax_radius = subhalos.vmax_radius

        spin = subhalos.spin
        velocity_dispersion = subhalos.velocity_dispersion

        if len(center_of_mass) == 0:
            return

        particles = snapshot.particles
        
        stars = particles.stars
        popIII = particles.sinks.popIII
        popII = particles.sinks.popII
        bhs = particles.sinks.black_holes
        gas = particles.gas

        catalog = GalaxyCatalog(snapshot=snap)

        for i in range(len(center_of_mass)):
            
            # Skip small proto-galactic systems
            if catalog_mass_type[i, 0] + catalog_mass_type[i, 4] + catalog_mass_type[i, 5] < 1e6:
                continue

            # Skip contaminated halos
            contamination = particles.contamination.inside(position[i], aperture_radius[i])

            if contamination.count > 0:
                continue

            stars_inside = stars.inside(position[i], aperture_radius[i])
            popIII_inside = popIII.inside(position[i], aperture_radius[i])
            popII_inside = popII.inside(position[i], aperture_radius[i])
            bhs_inside = bhs.inside(position[i], aperture_radius[i])
            gas_inside = gas.inside(position[i], aperture_radius[i])

            stellar_mass = stars_inside.total_mass + popIII_inside.total_mass + popII_inside.total_mass
            bh_mass = bhs_inside.total_mass

            primary_bh_id, primary_bh_mass = bhs_inside.largest_particle

            gas_cell_masses = gas_inside.masses
            gas_cell_densities = gas_inside.densities
            gas_cell_temperatures = gas_inside.temperatures
            gas_cell_positions = gas_inside.positions
            gas_cell_velocities = gas_inside.velocities
            gas_cell_metallicities = gas_inside.metallicities

            
            galaxy = Galaxy(
                
                snapshot = snap,
                fof_index = fof_indexes[i],
                subhalo_index = i,

                position = position[i],
                velocity = velocity[i],
                aperture_radius = aperture_radius[i],
                center_of_mass = center_of_mass[i],

                catalog_mass = catalog_mass[i],
                catalog_mass_type = catalog_mass_type[i, :],
                halfmass_radius = halfmass_radius[i],
                halfmass_radius_type = halfmass_radius_type[i, :],
                vmax = vmax[i],
                vmax_radius = vmax_radius[i],
                spin = spin[i],
                velocity_disperion = velocity_dispersion[i],

                gas_mass = gas_cell_masses.sum(),
                gas_cell_masses = gas_cell_masses,
                gas_cell_density = gas_cell_densities,
                gas_cell_temperature = gas_cell_temperatures,
                gas_cell_positions = gas_cell_positions,
                gas_cell_velocities = gas_cell_velocities,
                gas_cell_metallicity = gas_cell_metallicities,

                stellar_particle_types= np.concatenate([
                    np.full(len(stars_inside.ids), CS.POPII_OLD, dtype=np.uint8),
                    np.full(len(popIII_inside.ids), CS.POPIII, dtype=np.uint8),
                    np.full(len(popII_inside.ids), CS.POPII_YOUNG, dtype=np.uint8),
                ]),

                stellar_particle_ids = np.concatenate([
                    stars_inside.ids,
                    popIII_inside.ids,
                    popII_inside.ids,
                ]),

                stellar_particle_masses = np.concatenate([
                    stars_inside.masses,
                    popIII_inside.masses,
                    popII_inside.masses,
                ]),

                stellar_particle_positions = np.concatenate([
                    stars_inside.positions,
                    popIII_inside.positions,
                    popII_inside.positions,
                ]),

                stellar_particle_velocities = np.concatenate([
                    stars_inside.velocities,
                    popIII_inside.velocities,
                    popII_inside.velocities,
                ]),

                stellar_particle_lifetimes = np.concatenate([
                    stars_inside.stellar_lifetimes,
                    popIII_inside.stellar_lifetimes,
                    popII_inside.stellar_lifetimes,
                ]),

                stellar_particle_formation_times = np.concatenate([
                    stars_inside.formation_times,
                    popIII_inside.formation_times,
                    popII_inside.formation_times,
                ]),

                stellar_mass = stellar_mass,

                bh_mass = bh_mass,
                primary_bh_mass = primary_bh_mass,
                primary_bh_id = primary_bh_id,
                bh_particle_ids = bhs_inside.ids,
                bh_particle_masses = bhs_inside.masses,
                bh_separations = bhs_inside.distances_to(position[i]),
                bh_particle_positions = bhs_inside.positions,
                bh_particle_velocities = bhs_inside.velocities,
            )

            catalog.galaxies.append(galaxy)

        self.sim.writer.write_pickle(outfile, catalog.to_dict())

     


'''
class GalaxyCatalog:

    def __init__(self, sim):

        self.sim = sim
        self.base = sim.base
        self.binary_reader = sim.binary_reader     

    def filename(self, snap):
        return f"Galaxies_{snap:03d}.pkl"
    
    def build_catalog(self, snap):

        outfile = self.filename(snap)

        

        if os.path.exists(os.path.join(self.base, outfile)):
            print(f"Galaxy catalog already exists fro snap {snap}")
            return
        
        snapshot = self.sim.snapshot(snap)
        groups = snapshot.groups
        
        if groups.header["Ngroups_ThisFile"] == 0:
            print(f"There are no halos in the given snapshot. Skipping it.")
            return
        
        subhalos = groups.subhalos

        mass = subhalos["SubhaloMass"] * 1e10 / CS.hubble_parameter
        pos = subhalos["SubhaloCM"] 
        rad = subhalos["SubhaloHalfmassRad"]
        mtype = subhalos["SubhaloMassType"] * 1e10 / CS.hubble_parameter

        if len(mass) == 0:
            return
        

        YT = snapshot.yt
        
        YT.fields.metallicity()

        particles = snapshot.particles
        stars = particles.stars
        sinks = particles.sinks
        bhs = sinks.black_holes
        popIII = sinks.popIII
        popII = sinks.popII


        galaxies = {}

        for i in range(len(mass)):

            if mtype[i, 5] < 1e5: 
                continue

            contamination = particles.contamination.inside(pos[i], rad[i])
            if contamination.count > 0:
                continue

            stars_inside = stars.inside(pos[i], rad[i])
            bhs_inside = bhs.inside(pos[i], rad[i])
            popIII_inside = popIII.inside(pos[i], rad[i])
            popII_inside = popII.inside(pos[i], rad[i])

            bhmass = bhs_inside.total_mass
            bhmasses = bhs_inside.masses
            bhids = bhs_inside.ids
            mid, mmax = bhs_inside.largest_particle
        
            stellar = stars_inside.total_mass + popIII_inside.total_mass + popII_inside.total_mass
            sep = bhs_inside.distances_to(pos[i])

            print(bhmass, stellar)
            
            region = YT.sphere(pos[i], rad[i], units="code_length")

            gasmass, coldgasmass, hotgasmass = region.analysis.gas_mass
            
            Zavg = region.analysis.mean_metallicity

            print(gasmass, Zavg)

            galaxies[i] = dict(
                pos    = pos[i],
                radius = rad[i],
                mass   = mass[i],

                metallicity   = Zavg,
                stellar_mass  = stellar,
                gas_mass      = gasmass,
                cold_gas_mass = coldgasmass,
                hot_gas_mass  = hotgasmass,
                
                bh_mass         = bhmass,
                most_massive_bh = mmax,
                most_massive_id = mid,
                bh_ids          = bhids,
                bh_masses       = bhmasses,
                separation      = sep,
            )
        self.sim.writer.write_pickle(outfile, galaxies)
'''