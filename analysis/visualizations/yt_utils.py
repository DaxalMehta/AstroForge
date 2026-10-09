import yt
import numpy as np
import os
from scipy.spatial import cKDTree
import logging
logging.getLogger('yt').setLevel(logging.WARNING)

from ..core.utils import unwrap

class YTSession:

    def __init__(self, ds):
        
        self.ds = ds
        self._ad = None

        self.redshift = ds.current_redshift

        self.fields = FieldManager(self)
    
    @property
    def ad(self):
        if self._ad is None:
            self._ad = self.ds.all_data()
        return self._ad   

    def sphere(self, center, radius, velocity=None, units="kpc"):

        center = self.ds.arr(center, "code_length")
        radius = self.ds.quan(radius, units)
        if velocity is not None:
            velocity = self.ds.arr(velocity, "code_velocity")

        data = self.ds.sphere(center, radius)

        return Region(self, data, center=center, width=radius, velocity=velocity, geometry="sphere")

class FieldManager:

    def __init__(self, session):
        
        self.session = session
        self.ds = session.ds
    
    def _metallicity_field(self, field, data):
        Metallicity = np.divide(data[("PartType0", "SneTracerField")], data[("PartType0", "Masses")].v)
        Metallicity = np.clip(Metallicity, 1e-10, None)
        return Metallicity
    
    def metallicity(self):
        self.ds.add_field(("PartType0", "Metallicity"), units="", sampling_type="particle", function=self._metallicity_field)

    def _radial_velocity_field(self, field, data):
        center = self.center
        velocity = self.velocity.to("km/s")

        pos = data[("PartType0", "Coordinates")]
        vel = data[("PartType0", "Velocities")].to("km/s")

        rel_vel = vel - velocity 

        r = pos - center
        radius = np.linalg.norm(r, axis=1) * pos.units
        #radius = np.maximum(radius, 1e-30*radius.units)

        vr = np.sum(r * rel_vel, axis=1) / radius

        return vr

    def radial_velocity(self, center, velocity):
        self.center = center
        self.velocity = velocity

        self.ds.add_field(("PartType0", "RadialVelocity"), function=self._radial_velocity_field, units="km/s", sampling_type="particle")


    def dm_cic(self):
        self.ds.add_deposited_particle_field(("PartType1", "Masses"), method="cic", weight_field=None)

class Region:

    def __init__(self, session, data, center=None, width=None, velocity=None, geometry="generic"):
        
        self.session = session
        self.data = data
        
        self.center = center
        self.width = width
        self.velocity = velocity
        
        self.geometry = geometry
        
        self.analysis = RegionAnalysis(self)
    
    def plot(self):
            return Plot(self)

class RegionAnalysis:

    def __init__(self, region):
        
        self.region = region
        self.session = region.session
        self.data = region.data

        self.center = region.center
        self.width = region.width
      
    def contamination(self):

        if "PartType2" in self.session.ds.fields:
            Mass2 = self.data[("PartType2", "Masses")].in_units("Msun").v

            if np.sum(Mass2) > 0:
                print(f"Warning: Contamination particles found in the region Contanimation Mass = {np.sum(Mass2):.2e} Msun")
   
    @property
    def gas_mass(self):
        
        mass = self.data[("PartType0", "Masses")].in_units("Msun").value
        rho = self.data[("gas", "number_density")].value
        T = self.data[("gas", "temperature")].value

        mask = (T < 1e4) & (rho > 1) # Cold and Dense gas

        total_gas_mass = mass.sum()
        cold_gas_mass = mass[mask].sum()
        hot_gas_mass = total_gas_mass - cold_gas_mass

        return total_gas_mass, cold_gas_mass, hot_gas_mass

    @property
    def mean_metallicity(self):
        Z = self.data[("PartType0", "Metallicity")].v
        return np.mean(Z)

    @property
    def velocity_dispersion(self):

        mass = self.data[("PartType0", "Masses")].in_units("Msun").v
        vel = self.data[("PartType0", "Velocities")].in_units("km/s").v

        v_com = np.average(vel, axis=0, weights=mass)

        sigma = np.sqrt(np.sum(mass * np.sum((vel - v_com)**2, axis=1))/ np.sum(mass))

        return sigma


    def find_bh_and_star(self, id_binary, type_binary):
        id5 = self.data[("PartType5", "ParticleIDs")].value
        mass5 = self.data[("PartType5", "Masses")].in_units("Msun").value
        pos5 = self.data[("PartType5", "Coordinates")]
        dist = np.linalg.norm((pos5 - self.center), axis=1)

        sorter = np.argsort(id_binary)
        idx = np.searchsorted(id_binary, id5, sorter=sorter)
        mask = (idx >= 0) & (idx < len(type_binary))
        Type = np.zeros(len(id5))
        Type[mask] = type_binary[sorter[idx[mask]]]

        is_bh = Type == 3
        is_stellar = np.isin(Type, [0, 2])

        BHMass = np.sum(mass5[is_bh])
        StellarMass = np.sum(mass5[is_stellar])

        if "PartType4" in self.ds.fields:
            StellarMass += self.data[("PartType4", "Masses")].in_units("Msun").sum().to_value()
        
        ID = [int(unwrap(pid)) for pid in id5[is_bh]]
        HeavySeedMass = [unwrap(pmass) for pmass in mass5[is_bh]]
        Separation = [unwrap(pdist) for pdist in dist[is_bh]]
        
        if np.any(is_bh):
            i_max = np.argmax(mass5[is_bh])
            MostMassiveSeed = mass5[is_bh][i_max]
            MostMassiveID = id5[is_bh][i_max]
        else:
            MostMassiveSeed, MostMassiveID = 0.0, -1
        
        return BHMass, StellarMass, ID, HeavySeedMass, Separation, MostMassiveSeed, int(MostMassiveID)

class Plot:

    def __init__(self, region):
        
        self.region = region
        self.session = region.session
        
        @property
        def data(self):
            return self.region.data

        @property
        def center(self):
            return self.region.center

        @property
        def width(self):
            return self.region.width
        
        self.plot = None
        self.field = None
        self.weight_field = None
        self.normal = None
        self._frb = None

    
    @property
    def annotate(self):
        return Annotator(self)
  
    ##############################
    # Plot Constructors
    ##############################

    def _projection(self, field, weight_field=("gas", "number_density"), normal="x", north_vector=None):
        
        self.field, self.weight_field, self.normal = field, weight_field, normal
        self.plot = yt.ProjectionPlot(self.session.ds, normal, field, center=self.region.center, width=self.region.width, weight_field=weight_field, data_source=self.region.data)

        return self

    def _projection2(self, field, weight_field=("gas", "number_density"), normal="x", north_vector=None):

        self.field, self.weight_field, self.normal = field, weight_field, normal
        self.plot = yt.OffAxisProjectionPlot(self.session.ds, normal, field, center=self.region.center, width=self.region.width, weight_field=weight_field, data_source=self.region.data, north_vector=north_vector)

        return self

    def _particle_projection(self, field, weight_field=("PartType0", "Masses"), normal="x", deposition="cic", north_vector=None):

        self.field, self.weight_field, self.normal = field, weight_field, normal
        self.plot = yt.ParticleProjectionPlot(self.session.ds, normal, field, center=self.region.center, width=self.region.width, weight_field=weight_field, deposition=deposition, data_source=self.region.data, north_vector=north_vector)

        return self

    def _slice(self, field, weight_field=("gas", "number_density"), normal="x"):
        self.field, self.weight_field, self.normal = field, weight_field, normal

        self.plot = yt.SlicePlot(self.session.ds, normal, field, center=self.region.center, width=self.region.width, data_source=self.region.data)

        return self
    
    def phase(self):
        self.plot = yt.PhasePlot(self.data, ("gas", "number_density"), ("gas", "temperature"), ("gas", "mass"))
        self.field = ("gas", "number_density")
        self.weight_field = ("gas", "mass")
        self.normal = None

        return self
    
    def _particles(self, particle_type, cmap="grey", normal="x"):
        
        self.field = (particle_type, "Masses")
        self.weight_field = (particle_type, "Masses")
        self.normal = normal

        if normal == "x":
            ax1, ax2 = "y", "z"
        elif normal == "y":
            ax1, ax2 = "x", "z"
        elif normal == "z":
            ax1, ax2 = "x", "y"
        
        self.plot = yt.ParticlePlot(self.session.ds, (particle_type, "particle_position_"+ax1), (particle_type, "particle_position_"+ax2), (particle_type, "Masses"), center=self.region.center, width=self.region.width, data_source=self.region.data)

        self.plot.set_unit((particle_type, "Masses"), "Msun")
        self.plot.set_cmap((particle_type, "Masses"), cmap)
        self.plot.set_background_color((particle_type, "Masses"))

        return self
    
    
    ##############################
    #Convenience Wrappers
    ##############################

    def density(self, normal="x"):
        return self._projection(("gas", "number_density"), normal=normal)
    
    def temperature(self, normal="x", north_vector=None):
        return self._projection2(("gas", "temperature"), normal=normal, north_vector=north_vector)
    
    def metallicity(self, normal="x"):
        return self._projection(("PartType0", "Metallicity"), weight_field=("PartType0", "Density"), normal=normal)
    
    def radial_velocity(self, normal="x"):
        #return self._projection(("PartType0", "RadialVelocity"), weight_field=("PartType0", "Masses"), normal=normal)
        return self._slice(("PartType0", "RadialVelocity"), normal=normal)
        #return self._particle_projection(("PartType0", "RadialVelocity"))
    
    def dm(self, normal="x"):
        return self._particles("PartType1", normal=normal)
    
    def stars(self, normal="x"):
        return self._particles("PartType4", normal=normal)

    def stars2(self, camera_direction, north_vector=None, deposition="cic"):
        return self._particle_projection(field=("PartType4", "Masses"), weight_field=("PartType4", "Masses"), normal=camera_direction, deposition=deposition, north_vector=north_vector)


    def stars_camera(
        self,
        camera_direction,
        north_vector=(0.0, 0.0, 1.0),
        cmap="grey",
        marker_size=1.0,
    ):
        """
        Project stellar particles onto an arbitrary camera plane.

        Parameters
        ----------
        camera_direction : array-like
            Direction from the galaxy toward the camera.

        north_vector : array-like
            Direction that defines "up" in the resulting image.

        cmap : str
            Matplotlib colormap.

        marker_size : float
            Size of stellar particles.

        Returns
        -------
        Plot
            Plot object with a matplotlib figure.
        """

        camera_direction = np.asarray(camera_direction, dtype=float)
        north_vector = np.asarray(north_vector, dtype=float)

        # Normalize camera direction
        camera_direction /= np.linalg.norm(camera_direction)

        # Normalize north vector
        north_vector /= np.linalg.norm(north_vector)

        # Check that vectors are not parallel
        if np.isclose(abs(np.dot(camera_direction, north_vector)), 1.0):
            raise ValueError(
                "north_vector cannot be parallel to camera_direction."
            )

        # ------------------------------------------------------------
        # Construct camera basis
        # ------------------------------------------------------------

        right_vector = np.cross(north_vector, camera_direction)
        right_vector /= np.linalg.norm(right_vector)

        up_vector = np.cross(camera_direction, right_vector)
        up_vector /= np.linalg.norm(up_vector)

        # ------------------------------------------------------------
        # Get stellar particles
        # ------------------------------------------------------------

        positions = self.region.data[("PartType4", "Coordinates")]

        positions = positions.in_units("code_length").v

        center = self.region.center.in_units("code_length").v

        relative_position = positions - center

        # ------------------------------------------------------------
        # Project onto camera plane
        # ------------------------------------------------------------

        x = np.dot(relative_position, right_vector)
        y = np.dot(relative_position, up_vector)

        # ------------------------------------------------------------
        # Create matplotlib figure
        # ------------------------------------------------------------

        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(8, 8))

        ax.scatter(
            x,
            y,
            s=marker_size,
            c="white",
            cmap=cmap,
            linewidths=0,
        )

        # ------------------------------------------------------------
        # Set camera width
        # ------------------------------------------------------------

        width = self.region.width.in_units("code_length").v

        ax.set_xlim(-width / 2, width / 2)
        ax.set_ylim(-width / 2, width / 2)

        ax.set_aspect("equal")

        ax.set_axis_off()

        fig.savefig(f"Frame1")
        #self.plot = fig

        #self.field = ("PartType4", "Masses")
        #self.weight_field = ("PartType4", "Masses")
        #self.normal = camera_direction

        return self

    ##############################
    # Styling
    ############################## 
    
    def buff_size(self, buff_size):
        self.plot.set_buff_size(buff_size)
        return self
    
    def cmap(self, cmap):
        self.plot.set_cmap(self.field, cmap)
        return self

    def background(self):
        self.plot.set_background_color(self.field)
        return self
    
    def units(self, units):
        self.plot.set_unit(self.field, units)
        return
    
    def limits(self, lowlim, highlim):
        self.plot.set_zlim(self.field, lowlim, highlim)
        return self
    
    def beauty(self):
        self.plot.hide_colorbar()
        self.plot.hide_axes(draw_frame=True)
        return self
    
    ##############################
    # Saving
    ##############################

    def save(self, base="", prefix=""):
        
        if self.plot is None or self.field is None:
            raise RuntimeError("Create a plot before saving it.")
        if base:
            base = "Images/" + base
            os.makedirs(base, exist_ok=True)
        fname = f"{prefix}_{self.field[1].replace(' ', '_')}.png"
        path = os.path.join(base, fname)
        self.plot.save(path)
        print(f"Plot created at {path}")
        return path

    ##############################
    # frb
    ##############################
    @property
    def frb(self):

        if self._frb is None:

            self._frb = self.plot.data_source.to_frb(
                self.region.width,
                512
            )

        return self._frb
    
    def array(self):
        return np.asarray(self.frb[self.field])

class Annotator:

    def __init__(self, Plot):
        
        self.region = Plot.region
        self.session = Plot.session

        self.plot = Plot.plot
        self.center = Plot.region.center
        self.radius = Plot.region.width
    
    def scale(self):
        self.plot.annotate_scale(size_bar_args={'color': 'w'}, coeff=100, unit='pc', text_args={'size':15})
    
    def particles(self, positions, color, marker, alpha, s=1):

        for position in positions:
            self.plot.annotate_marker(position, coord_system="data", color=color, s=s, marker=marker, alpha=alpha)
        
    def old_stars(self, color):

        if "PartType4" in self.ds.fields:
            Pos4 = self.region[("PartType4", "Coordinates")]
            for position in Pos4:
                self.plot.annotate_marker(position, coord_system="data", color=color, s=1, marker="o", alpha=0.3)
    
    def blackholes(self, positions, masses, color):

        for position, mass in zip(positions, masses):

            if mass < 1e3:
                self.plot.annotate_marker(position, coord_system="data", color=color, s=3, marker="o")
            elif mass > 1e3 and mass < 1e5:
                self.plot.annotate_marker(position, coord_system="data", color=color, s=4, marker="^")
            elif mass > 1e5:
                self.plot.annotate_marker(position, coord_system="data", color=color, s=5, marker="d")
    
    def halo(self, center, radius, units="kpc"):
        center = self.session.ds.arr(center, "code_length")
        radius = self.session.ds.quan(radius, units)
        
        self.plot.annotate_marker(center, coord_system="data", color="red", s=50, marker="x")
        self.plot.annotate_sphere(center, radius=radius, circle_args={"color": "red"})
    
    def text(self, locx, locy, text, fontsize=50, fontcolor="white"):
        self.plot.annotate_text((locx, locy), text, coord_system="axis", text_args={"color":fontcolor, "fontsize":fontsize, "fontweight":"bold"})
    