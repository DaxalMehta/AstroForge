import yt
import numpy as np
import os
from scipy.spatial import cKDTree
import logging
logging.getLogger('yt').setLevel(logging.WARNING)

from ..core.utils import unwrap

class YT:
    def __init__(self, ds, ad=None):
        self.ds = ds
        self.ad = ad if ad is not None else ds.all_data()
        self.region_tools = YTRegion(ds)
        self.plot = Plot(ds)

    def _plotter_for_current_region(self, plotter_class, normal="x", frb_output=False):
        """Create a plotter that uses the region currently selected on this YT object."""
        plotter = plotter_class(self.ds, normal=normal, frb_output=frb_output)
        if self.region_tools.region is not None:
            plotter.region = self.region_tools.region
            plotter.center = self.region_tools.center
            plotter.width = self.region_tools.width
        return plotter

    def new_plot(self, normal="x"):
        """Return an independent general-purpose plotter for the selected region."""
        return self._plotter_for_current_region(Plot, normal)

    def number_density_projection(self, normal="x", frb_output=False):
        """Return an independent number-density projection for the selected region."""
        return self._plotter_for_current_region(NumberDensityProjection, normal, frb_output)
    
    def temperature_projection(self, normal="x", frb_output=False):
        return self._plotter_for_current_region(TemperatureProjection, normal, frb_output)
    
    def metallicity_projection(self, normal="x", frb_output=False):
        return self._plotter_for_current_region(MetallicityProjection, normal, frb_output)
    
    def stellar_projection(self, normal="x", frb_output=False):
        return self._plotter_for_current_region(StellarProjectionPlot, normal, frb_output)
    
    def region(self, center, radius, units="kpc"):
        region = self.region_tools.SelectRegion(center, radius, units)
        self.plot.region = region
        self.plot.center = self.region_tools.center
        self.plot.width = self.region_tools.width
        return region  
    
    def calculate_high_res_volume(self):

        # -----------------------------
        # METHOD 1: Bounding box of high-res particles (PartType1)
        # -----------------------------

        pos1 = self.ad[("PartType1", "Coordinates")].to_value()

        xmin, ymin, zmin = np.min(pos1, axis=0)
        xmax, ymax, zmax = np.max(pos1, axis=0)

        V1 = (xmax - xmin) * (ymax - ymin) * (zmax - zmin)

        # -----------------------------
        # METHOD 2: Region not contaminated by PartType2
        # -----------------------------

        pos2 = self.ad[("PartType2", "Coordinates")].to_value()

        if len(pos2) == 0:
            # No contamination → full volume is high-res
            return V1, V1
        
        center = 0.5 * (np.max(pos1, axis=0) + np.min(pos1, axis=0))

        print(f"Center of the box is {center}")

        pos2_shift = pos2 - center
        x, y, z = pos2_shift[:, 0], pos2_shift[:, 1], pos2_shift[:, 2]

        def find_limits(coord):
            pos = coord[coord > 0]
            neg = coord[coord < 0]

            if len(pos) == 0 or len(neg) == 0:
                return None, None

            return np.min(pos), np.max(neg)

        xmax2, xmin2 = find_limits(x)
        ymax2, ymin2 = find_limits(y)
        zmax2, zmin2 = find_limits(z)

        # Handle edge cases
        if None in (xmax2, xmin2, ymax2, ymin2, zmax2, zmin2):
            V2 = 0.0
        else:
            V2 = (xmax2 - xmin2) * (ymax2 - ymin2) * (zmax2 - zmin2)

        tree = cKDTree(pos2)
        dist, _ = tree.query(pos1, k=1)
        r_max = np.min(dist)

        V3 = (4/3) * np.pi * r_max**3 

        return V1, V2, V3


class YTRegion:
    def __init__(self, ds):
        self.ds = ds
        self.fields = []
        self.region = None
        self.center = ds.domain_center
        self.width = "1000"
    
    def SelectRegion(self, center, radius, units="kpc"):
        if units == "code_length":
            radius =  radius / (1 + self.ds.current_redshift) / self.ds.hubble_constant
            units = "kpc"
            if radius > 1500:
                radius /= 1000
                units = "Mpc"
        self.center = self.ds.arr(center, "code_length")
        self.width = self.ds.quan(radius, units)
        self.region = self.ds.sphere(self.center, self.width)

        return self.region
    
    def _Metallicity(self, field, data):
        Metallicity = np.divide(data[("PartType0", "SneTracerField")], data[("PartType0", "Masses")].v)
        Metallicity = np.clip(Metallicity, 1e-10, None)
        return Metallicity
    
    def add_metallicity(self, weight_field=("PartType0", "Masses")):
        self.ds.add_field(("PartType0", "Metallicity"), units="", sampling_type="particle", function=self._Metallicity)
        field = ("PartType0", "Metallicity")
        self.fields.append((field, weight_field))

    def DMfields(self):
        self.ds.add_deposited_particle_field(("PartType1", "Masses"), method="cic", weight_field=None)

    def FindBHandStellarProperty(self, id_binary, type_binary):
        id5 = self.region[("PartType5", "ParticleIDs")].value
        mass5 = self.region[("PartType5", "Masses")].in_units("Msun").value
        pos5 = self.region[("PartType5", "Coordinates")]
        dist = np.linalg.norm((pos5 - self.region.center), axis=1)

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
            StellarMass += self.region[("PartType4", "Masses")].in_units("Msun").sum().to_value()
        
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
    
    def CheckContamination(self):
        if "PartType2" in self.ds.fields:
            Mass2 = self.region[("PartType2", "Masses")].in_units("Msun").value
            if np.sum(Mass2) > 0:
                print("Warning: Contamination particles found in the region. Contanimation Mass = %.2e Msun" % np.sum(Mass2))

    def FindGasMass(self):
        if self.region is None:
            raise RuntimeError("Select a region before calculating gas masses.")
        GasMass = self.region[("PartType0", "Masses")].in_units("Msun").value
        GasDensity = self.region[("gas", "number_density")].value
        GasTemp = self.region[("gas", "temperature")].value
        mask = (GasTemp < 1e4) & (GasDensity > 1) ## Cold and Dense gas
        TotalGasMass = GasMass.sum()
        ColdGasMass = GasMass[mask].sum()
        HotGasMass = TotalGasMass - ColdGasMass
        return TotalGasMass, ColdGasMass, HotGasMass

    def find_gas_mass(self):
        """PEP 8 alias for :meth:`FindGasMass`."""
        return self.FindGasMass()

    def get_frb_data(self, normal="x", buff_size=512):
        """TODO: Write a get_frb_data"""

class YTPlotter(YTRegion):
    def __init__(self, ds, normal="x"):
        super().__init__(ds)
        self.normal = normal
        self.plot = None
        self.field = None
        self.region = None
        self.center = ds.domain_center
        self.width = "1000"
    
    def reset(self):
        """Clear plot-specific state while preserving the selected region."""
        self.plot = None
        self.field = None
        self.fields = []
        return self

    def add_field(self, field, weight_field=("gas", "number_density")):
        self.fields = [(field, weight_field)]
        return self
    
    def annotate_particles(self, positions, color, marker, alpha):

        for position in positions:
            self.plot.annotate_marker(position, coord_system="data", color=color, s=1, marker=marker, alpha=alpha)
    
    def annotate_old_stars(self, color):
        
        if "PartType4" in self.ds.fields:
            Pos4 = self.region[("PartType4", "Coordinates")]
            for position in Pos4:
                self.plot.annotate_marker(position, coord_system="data", color=color, s=1, marker="o", alpha=0.3)
    
    def annotate_blackholes(self, positions, masses, color):

        for position, mass in zip(positions, masses):

            if mass < 1e3:
                self.plot.annotate_marker(position, coord_system="data", color=color, s=3, marker="o")
            elif mass > 1e3 and mass < 1e5:
                self.plot.annotate_marker(position, coord_system="data", color=color, s=4, marker="^")
            elif mass > 1e5:
                self.plot.annotate_marker(position, coord_system="data", color=color, s=5, marker="d")
    
    def plot_sink_particles(self, id_binary, type_binary, heavy_only=False, stars=False):
           
        id5 = self.region[("PartType5", "ParticleIDs")].value
        pos5 = self.region[("PartType5", "Coordinates")]
        mass5 = self.region[("PartType5", "Masses")].in_units("Msun").value
        sorter = np.argsort(id_binary)
        idx = np.searchsorted(id_binary, id5, sorter=sorter)

        Type = np.zeros(len(id5))
        mask = (idx < len(id_binary))
        Type[mask] = type_binary[sorter[idx[mask]]]
        color_map = {0: "cyan", 2: "gold", 3: "red"}
        size_map = {0: 2, 2: 2, 3: 5}
        for j in range(len(id5)):
            color = color_map.get(Type[j], None)
            size = size_map.get(Type[j], 1)
            if (color == "red" and mass5[j] > 1e6): color = "black"
            if heavy_only and color != "black": continue
            if color:
                self.plot.annotate_marker(pos5[j], coord_system="data", color=color, s=size, marker="o")
        return Type

    def plot_CMBH(self, Pos, Mass):

        self.plot.annotate_marker(Pos, coord_system="data", color="red", s=50, marker="d")
        #self.plot.annotate_text((0.05, 0.95), f"CMBH Mass: {unwrap(Mass):.2e} Msolar", coord_system="axis", text_args={"color":"white", "fontsize":50, "fontweight":"bold"})

    def plot_halo(self, halo_center, halo_radius, units="kpc"):
        if units == "code_length":
            halo_radius = halo_radius / (1 + self.ds.current_redshift) / self.ds.hubble_constant
            units = "kpc"
        radius = self.ds.quan(halo_radius, units)
        center = self.ds.arr(halo_center, "code_length")
        self.plot.annotate_marker(center, coord_system="data", color="red", s=50, marker="x")
        self.plot.annotate_sphere(center, radius=radius, circle_args={"color": "red"})
    
    def customize(self, cmap=None, background=False, lowlim=None, highlim=None, units=None, buff_size=None, beauty=False, scale=False):
        if self.plot is None or self.field is None:
            raise RuntimeError("Create a plot before customizing it.")
        plot = self.plot
        field = self.field
        if buff_size: plot.set_buff_size(buff_size)         
        if cmap: plot.set_cmap(field, cmap)
        if background == True: plot.set_background_color(field)
        if units: plot.set_unit(field, units)
        if lowlim is not None and highlim is not None: 
            if highlim <= lowlim:
                raise ValueError(f"Invalid limits: lowlim={lowlim}, highlim={highlim}")
            plot.set_zlim(field, lowlim, highlim)
            print(f"field= {field}, lowlim={lowlim}, highlim={highlim}")
        if scale == True:
            plot.annotate_scale(size_bar_args={'color': 'w'}, coeff=100, unit='pc', text_args={'size':15})  
        if beauty == True:
            plot.hide_colorbar()
            plot.hide_axes(draw_frame=True)           
                
    def annotate_text(self, locx, locy, text, fontsize=50, fontcolor="white"):
        self.plot.annotate_text((locx, locy), text, coord_system="axis", text_args={"color":fontcolor, "fontsize":fontsize, "fontweight":"bold"})
    
    def save_plots(self, base="", prefix=""):
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

class Plot(YTPlotter):
    def __init__(self, ds, normal="x", frb_output=False):
        super().__init__(ds, normal)

        self.frb_output = frb_output

    def ProjectionPlot(self):
        if not self.fields:
            raise RuntimeError("Add a field before creating a projection plot.")
        field, weight = self.fields[0]
        self.plot = yt.ProjectionPlot(self.ds, self.normal, field, center=self.center, width=self.width, weight_field=weight, data_source=self.region)
        self.field = field

        if self.frb_output == True:
            frb = self.plot.to_frb(self.width, 512)
            array = np.array(frb[field])           
            return array
        else:
            return None

    def SlicePlot(self):
        if not self.fields:
            raise RuntimeError("Add a field before creating a slice plot.")
        field, weight = self.fields[0]
        self.plot = yt.SlicePlot(self.ds, self.normal, field, center=self.center, width=self.width, weight_field=weight, data_source=self.region)
        self.field = field

        return self.plot

    def PhasePlot(self):
        self.plot = yt.PhasePlot(self.region, ("gas", "number_density"), ("gas", "temperature"), ("gas", "mass"))
        self.field = ("gas", "number_density")
        return self.plot
    
    def ParticlePlot(self, cmap="magma"):
        p = yt.ParticlePlot(self.ds, ("PartType1", "particle_position_y"), ("PartType1", "particle_position_z"), ("PartType1", "Masses"), center=self.center, width = self.width, data_source=self.region)

        p.set_unit(("PartType1", "Masses"), "Msun")
        p.set_cmap(("PartType1", "Masses"), cmap)
        p.set_background_color(("PartType1", "Masses"))

        self.plot = p
        self.field = ("PartType1", "particle_position_x")

        return self.plot

    def ParticleProjectionPlot(self, cmap="magma"):
        if not self.fields:
            raise RuntimeError("Add a field before creating a particle projection plot.")
        field, weight = self.fields[0]
        self.plot = yt.ParticleProjectionPlot(self.ds, self.normal, field, weight_field=weight, deposition="cic", center=self.center, width=self.width)
        self.field = field

        return self.plot

    def StellarParticlePlot(self, cmap="magma"):
        p = yt.ParticlePlot(self.ds, ("PartType4", "particle_position_y"), ("PartType4", "particle_position_z"), ("PartType4", "Masses"), center=self.center, width = self.width, data_source=self.region)

        p.set_unit(("PartType4", "Masses"), "Msun")
        p.set_cmap(("PartType4", "Masses"), cmap)
        p.set_background_color(("PartType4", "Masses"))

        self.field = ("PartType4", "particle_position_x")
        self.plot = p

        return self.plot

class NumberDensityProjection(Plot):
    """Convenience plotter for gas number-density projections.

    Instances are normally created with ``YT.number_density_projection()`` so
    they automatically inherit the region selected by ``YT.region()``.
    """

    FIELD = ("gas", "number_density")

    def __init__(self, ds, normal="x", frb_output=False):
        super().__init__(ds, normal=normal, frb_output=frb_output)
        self.add_field(self.FIELD, weight_field=self.FIELD)

        self.frb_output = frb_output
    
    def reset(self):
        """Reset the rendered plot and restore the number-density field."""
        super().reset()
        self.add_field(self.FIELD, weight_field=self.FIELD)
        return self

    def create_plot(self, cmap="viridis", lowlim=1e-1, highlim=1e5, background=True, beauty=True, scale=True, units=None, buff_size=None):
        """Create and style the projection, returning the yt plot object."""
        
        frb = self.ProjectionPlot()
        self.customize(cmap=cmap, lowlim=lowlim, highlim=highlim, background=background, beauty=beauty, scale=scale, units=units, buff_size=buff_size)
        if self.frb_output == True:
            return frb
        else:
            return None

class MetallicityProjection(Plot):
    """
    Creates a projection plot for the metallicity field.
    """

    FIELD = ("PartType0", "Metallicity")

    def __init__(self, ds, normal="x", frb_output=False):
        super().__init__(ds, normal, frb_output=frb_output)
        self.add_field(self.FIELD)

        self.frb_output = frb_output

    
    def reset(self):
        """Reset the rendered plot"""
        super().reset()
        self.add_field(self.FIELD)
        return self
    
    def create_plot(self, cmap="cividis", lowlim=1e-5, highlim=1e-2, background=True, beauty=True, scale=True, units=None, buff_size=None):
        """Create and style the projection, returning the yt plot object."""
        
        frb = self.ProjectionPlot()
        self.customize(cmap=cmap, lowlim=lowlim, highlim=highlim, background=background, beauty=beauty, scale=scale, units=units, buff_size=buff_size)
        
        if self.frb_output == True:
            return frb
        else:
            return None

class TemperatureProjection(Plot):
    """
    Creates a projection plot for the temperature field.
    """

    FIELD = ("gas", "temperature")

    def __init__(self, ds, normal="x", frb_output=False):
        super().__init__(ds, normal, frb_output=frb_output)
        self.add_field(self.FIELD)

        self.frb_output = frb_output
    
    def reset(self):
        """Reset the rendered plot"""
        super().reset()
        self.add_field(self.FIELD)
        return self
    
    def create_plot(self, cmap="magma", lowlim=1e2, highlim=5e5, background=True, beauty=True, scale=True, units=None, buff_size=None):
        """Create and style the projection, returning the yt plot object."""
        
        frb = self.ProjectionPlot()
        self.customize(cmap=cmap, lowlim=lowlim, highlim=highlim, background=background, beauty=beauty, scale=scale, units=units, buff_size=buff_size)

        if self.frb_output == True:
            return frb
        else:
            return None

class StellarProjectionPlot(Plot):
    """
    Creates a projection plot for the stellar particles field.
    """

    def __init__(self, ds, normal="x", frb_output=False):
        super().__init__(ds, normal, frb_output=frb_output)
    
    def reset(self):
        """Reset the rendered plot"""
        super().reset()
        return self
    
    def create_plot(self, cmap="grey", beauty=True, scale=True, units=None, buff_size=512):
        """Create and style the projection, returning the yt plot object."""
        
        self.StellarParticlePlot(cmap=cmap)
        self.customize(beauty=beauty, scale=scale, units=units, buff_size=buff_size)
        return self.plot