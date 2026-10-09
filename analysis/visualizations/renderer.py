import numpy as np
from scipy.ndimage import gaussian_filter
from matplotlib import colors
import matplotlib.pyplot as plt

from ..core.constants import CS
from ..core.utils import find_cosmic_time


class Scene:
    def __init__(self, sim, snap, center, radius, buff_size, projection_axis="x", camera_direction=None, north_vector=None):

        """
        Create a scene object.

        Parameters
        ----------
        sim : Simulation object
            Simulation on which we create the scene
        snap: Snapshot number
            The snapshot on which we create the scene
        center: array of length 3
            The center of the scene
        radius: float
            The extent of the scene
        buff_size: int
            The length of the pixel grid
        projection_axis: 'x', 'y', or 'z'
            The normal direction of the projection plane
        camera_direction: [float, float, float]
            The direction the camera is pointing at
        north_vector: [float, float, float]
            The north direction for the scene
        
        """
        
        self.sim = sim
        self.buff_size = buff_size

        self.projection_axis = projection_axis
        self.camera_direction = camera_direction
        self.north_vector = north_vector

        self.snap = snap
        self.center = center
        self.radius = radius
        
        self.snapshot = sim.snapshot(snap)
        self.particles = self.snapshot.particles
        self.sink_snap = self.snapshot.sink_snap
        self.snapshot_time = find_cosmic_time(self.sink_snap.time)

        self.render_tools = RenderTools(self)

    def update_camera(self, camera_direction=None, north_vector=None):
        if camera_direction is not None:
            self.camera_direction = np.asarray(camera_direction, dtype=float)
        if north_vector is not None:
            self.north_vector = np.asarray(north_vector, dtype=float)

        self.render_tools.setup_camera()

    def stellar_renderer(self):
        return StellarRenderer(self)
    
    def bh_renderer(self):
        return BHRenderer(self)

    def gas_renderer(self):
        return GasRenderer(self)

        
class RenderTools:
    def __init__(self, Scene):
        
        self.scene = Scene
        #self._setup_axes(Scene.projection_axis)
        self.setup_camera()

    def _setup_axes(self, axis):
        
        axis = axis.lower()
        if axis == "x":
            self.x_ax, self.y_ax, self.los_ax = 1, 2, 0
        elif axis == "y":
            self.x_ax, self.y_ax, self.los_ax = 0, 2, 1
        elif axis == "z":
            self.x_ax, self.y_ax, self.los_ax = 0, 1, 2
        else:
            raise ValueError("projection_axis must be 'x', 'y', or 'z'")

    def setup_camera(self):

        if self.scene.camera_direction is None:

            axis = self.scene.projection_axis.lower()

            camera_map = {
                "x": np.array([1.0, 0.0, 0.0]),
                "y": np.array([0.0, 1.0, 0.0]),
                "z": np.array([0.0, 0.0, 1.0]),
            }

            north_map = {
                "x": np.array([0.0, 0.0, 1.0]),
                "y": np.array([0.0, 0.0, 1.0]),
                "z": np.array([0.0, 1.0, 0.0]),
            }

            camera = camera_map[axis]
            north = north_map[axis]

        else:
            camera = np.asarray(self.scene.camera_direction, dtype=float)
            north = np.asarray(self.scene.north_vector, dtype=float)
        
        # Normalize
        camera /= np.linalg.norm(camera)
        
        #Remove the component of north paralle to camera.
        north -= np.dot(north, camera) * camera
        north /= np.linalg.norm(north)
        
        # They cannot be parallel
        if np.isclose(abs(np.dot(camera, north)), 1.0):
            raise ValueError("north_vector cannot be parallel to camera_direction")
        
        # Horizontal direction in image
        right = np.cross(north, camera)
        right /= np.linalg.norm(right)
        
        # Vertical direction in image
        up = np.cross(camera, right)
        up /= np.linalg.norm(up)
        
        self.camera_direction = camera
        self.north_vector = north
        self.right_vector = right
        self.up_vector = up
        
    def _coordinate_to_pixel(self, pos):
        """Convert coordinates in projection planes to pixel"""

        pos = np.asarray(pos)
        rel_pos = pos - self.scene.center
        
        #dx = rel_pos[self.x_ax]
        #dy = rel_pos[self.y_ax]

        #dx = np.dot(rel_pos, self.right_vector)
        #dy = np.dot(rel_pos, self.up_vector)
        dx = rel_pos @ self.right_vector
        dy = rel_pos @ self.up_vector

        xpix = ((dx + self.scene.radius) / (2 * self.scene.radius) * (self.scene.buff_size - 1)).astype(int)
        ypix = ((dy + self.scene.radius) / (2 * self.scene.radius) * (self.scene.buff_size - 1)).astype(int)      

        return xpix, ypix

    def _post_processing(self, R, G, B):
        """Normalize and apply asinh stretch"""
        RGB = np.dstack((R, G, B))

        nonzero = RGB[RGB>0]
        if nonzero.size == 0:
            print("There is no light in the projection. No point in plotting it.")
            return RGB
        
        stretch = np.percentile(nonzero, 98)
        if stretch <= 0:
            stretch = nonzero.max() if nonzero.max() > 0 else 1.0
        
        RGB = np.arcsinh(RGB/stretch)
        RGB /= RGB.max() if RGB.max() > 0 else 1.0

        return RGB

    def _pixel_radius(self, r):
        return int(self.scene.buff_size * r / (self.scene.radius * self.scene.sink_snap.time * 1e3 / CS.hubble_parameter))

    @staticmethod
    def _add_gaussian(R, G, B, row, col, L, color, sigma_core, sigma_halo, halo_fraction, N):
        ''' Add core + halo Gaussians to matrices R, G, B at (row, col)'''

        # Define a small patch around the center for efficiency
        # Use 4 * sigma to capture 99.99% of the flux

        rad_pix = int(np.ceil(4.0 * max(sigma_core, sigma_halo)))
        r0 = max(row - rad_pix, 0)
        r1 = min(row + rad_pix + 1, N)
        c0 = max(col - rad_pix, 0)
        c1 = min(col + rad_pix + 1, N)

        if r0 >= r1 or c0 >= c1:
            return
        
        # Coodinate grids
        yy, xx = np.mgrid[r0:r1, c0:c1]
        dy = yy - row
        dx = xx - col

        r2 = dx**2 + dy**2

        #Flux split
        L_halo = halo_fraction * L
        L_core = (1.0 - halo_fraction) * L

        # 2D Gaussian normalization sum(G) = 1

        core = np.exp(-0.5 * r2 / (sigma_core**2))
        core_sum = core.sum() if core.size > 0 else 1.0
        core *= (L_core / core_sum)

        halo = np.exp(-0.5 * r2 / (sigma_halo**2))
        halo_sum = halo.sum() if halo.size > 0 else 1.0
        halo *= (L_halo / halo_sum)

        Gtot = core + halo

        # Add to the color planes

        R[r0:r1, c0:c1] += Gtot * color[0]
        G[r0:r1, c0:c1] += Gtot * color[1]
        B[r0:r1, c0:c1] += Gtot * color[2]  

class StellarRenderer:
    
    def __init__(self, Scene):
        
        self.scene = Scene
        self.render_tools = Scene.render_tools
        self.sim = Scene.sim
        self.particles = Scene.particles

        stars = self.particles.stars.inside(self.scene.center, self.scene.radius)
        popii = self.particles.sinks.popII.inside(self.scene.center, self.scene.radius)

        self.stars = stars
        self.popii = popii

        # Calculate properties
        self.star_formation_time = find_cosmic_time(np.asarray(stars.formation_times))
        self.popii_formation_time = find_cosmic_time(np.asarray(popii.formation_times))

        self.star_age = (self.scene.snapshot_time - self.star_formation_time) / 1e6
        self.popii_age = (self.scene.snapshot_time - self.popii_formation_time) / 1e6

        # Color and Luminosity
        self.star_color, self.star_luminosity = self.color_by_age(self.star_age)
        self.popii_color, self.popii_luminosity = self.color_by_age(self.popii_age)


    def color_by_age(self, age):
        """
        Returns:
            color: RGB values
            luminosity: scaled value
        """

        age_grid = np.array([0, 10, 30, 75, 150, 367], dtype=float)

        rgb_grid = np.array([
            [0.15, 0.35, 1.00],   # newborn blue
            [0.35, 0.65, 1.00],   # blue
            [0.75, 0.90, 1.00],   # cyan
            [1.00, 0.95, 0.70],   # pale yellow
            [1.00, 0.60, 0.20],   # orange
            [1.00, 0.15, 0.10],   # red
        ])

        age = np.asarray(age)

        #age = np.clip(age, age_grid[0], age_grid[-1])

        r = np.interp(age, age_grid, rgb_grid[:, 0])
        g = np.interp(age, age_grid, rgb_grid[:, 1])
        b = np.interp(age, age_grid, rgb_grid[:, 2])

        color = np.column_stack((r, g, b))

        luminosity = 1/(1 + (age/60)**2)#np.exp(-age/50) #Young stars are brighter

        return color, luminosity
    
    def render(self):

        R = np.zeros((self.scene.buff_size, self.scene.buff_size), dtype=float)
        G = np.zeros_like(R)
        B = np.zeros_like(R)

        #stars = self.particles.stars
        #sinks = self.particles.sinks
        #popii = sinks.popII

        #stars = stars.inside(self.scene.center, self.scene.radius)
        #popii = popii.inside(self.scene.center, self.scene.radius)

        xpix, ypix = self.render_tools._coordinate_to_pixel(self.stars.positions)
        #formation_times = np.asarray(stars.formation_times)
        #formation_time = find_cosmic_time(formation_times)
        #age = (self.scene.snapshot_time - formation_time) / 1e6
        #color, luminosity = self.color_by_age(age)
        luminosity = self.star_luminosity * self.stars.masses
        
        np.add.at(R, (ypix, xpix), luminosity * self.star_color[:, 0])
        np.add.at(G, (ypix, xpix), luminosity * self.star_color[:, 1])
        np.add.at(B, (ypix, xpix), luminosity * self.star_color[:, 2])

        xpix, ypix = self.render_tools._coordinate_to_pixel(self.popii.positions)
        #formation_times = np.asarray(popii.formation_times)
        #formation_time = find_cosmic_time(formation_times)
        #age = (self.scene.snapshot_time - formation_time) / 1e6
        #color, luminosity = self.color_by_age(age)
        luminosity = self.popii_luminosity * self.popii.masses
        np.add.at(R, (ypix, xpix), luminosity * self.popii_color[:, 0])
        np.add.at(G, (ypix, xpix), luminosity * self.popii_color[:, 1])
        np.add.at(B, (ypix, xpix), luminosity * self.popii_color[:, 2])

        '''
        N = stars.count
        for i in range(N):
            p4 = stars.positions[i]
            m4 = stars.masses[i]

            xpix, ypix = self.render_tools._coordinate_to_pixel(p4)

            formation_time = find_cosmic_time(stars.formation_times[i])
            age = (self.scene.snapshot_time - formation_time) / 1e6 # in Myr

            color, lum = self.color_by_age(age)

            if not np.all(np.isfinite(color)):
                print("Bad color:", color, age)
                continue

            if not np.isfinite(lum):
                print("Bad luminosity:", lum, age)
                continue

            L = lum * m4

            R[ypix, xpix] += L * color[0]
            G[ypix, xpix] += L * color[1]
            B[ypix, xpix] += L * color[2]
        
        
        N = popii.count
        for i in range(N):
            p4 = popii.positions[i]
            m4 = popii.masses[i]

            xpix, ypix = self.render_tools._coordinate_to_pixel(p4)

            formation_time = find_cosmic_time(popii.formation_times[i])
            age = (self.scene.snapshot_time - formation_time) / 1e6 # in Myr

            color, lum = self.color_by_age(age)

            L = lum * m4

            R[ypix, xpix] += L * color[0]
            G[ypix, xpix] += L * color[1]
            B[ypix, xpix] += L * color[2]
        '''
        
        R = gaussian_filter(R, 0.4)
        G = gaussian_filter(G, 0.4)
        B = gaussian_filter(B, 0.4)
        RGB = self.render_tools._post_processing(R, G, B)

        return RGB
                 
class BHRenderer:
    def __init__(self, Scene, bh_color=np.array([0.9, 0.95, 1.0])):
        
        self.scene = Scene
        self.sim = Scene.sim
        self.render_tools = Scene.render_tools

        self.black_holes = self.scene.particles.sinks.black_holes
        self.black_holes = self.black_holes.inside(self.scene.center, self.scene.radius)
        self.black_hole_histories = self.sim.history.black_holes

        self.bh_color = bh_color

        self.luminosities = self.find_bh_luminosity()
        #self.positions = self.find_bh_pos()
    
    def find_bh_luminosity(self):
        
        luminosity = []
        bh_ids = self.black_holes.ids
        for bid in bh_ids:
            if bid not in self.black_hole_histories: 
                luminosity.append(0)
                continue
            
            history = self.black_hole_histories[bid]
            times = history.dump_times
            masses = np.array(history.masses) * 1e10 / CS.hubble_parameter
            
            if len(masses) < 2:
                luminosity.append(0)
                continue

            rate = (masses[-1] - masses[-2]) / (find_cosmic_time(times[-1]) - find_cosmic_time(times[-2])) * CS.Msun_to_g / CS.yr_to_s

            luminosity.append(0.1 * rate * CS.c_cgs**2)
        
        return np.array(luminosity)
    
    def find_bh_pos(self):

        #if self.black_holes.positions:
        #    return self.black_holes.positions
        pos = []
        for bid in self.scene.bh_ids:
            if bid not in self.sink_particles: 
                pos.append([np.nan, np.nan, np.nan])
                continue

            evo = self.sink_particles[bid]["evolution"]
            times = np.array(list(sorted(evo.keys())))

            pos.append(evo[np.max(times)]["Pos"][0])
        
        return np.array(pos)

    def render(self):

        R = np.zeros((self.scene.buff_size, self.scene.buff_size), dtype=float)
        G = np.zeros_like(R)
        B = np.zeros_like(R)
        
        for lum, pos in zip(self.luminosities, self.black_holes.positions):
            
            xpix, ypix = self.render_tools._coordinate_to_pixel(pos)

            brightness = lum + 1e30

            self.render_tools._add_gaussian(R, G, B, ypix, xpix, brightness, self.bh_color, 1, 3, 0.3, self.scene.buff_size)
        
        
        RGB = self.render_tools._post_processing(R, G, B)

        return RGB

class GasRenderer:

    def __init__(self, Scene):
    
        self.scene = Scene
        self.sim = Scene.sim
        self.render_tools = Scene.render_tools

        gas = self.scene.particles.gas
        self.gas = gas.inside(Scene.center, Scene.radius)
        count = self.gas.count
        self.xpix = []
        self.ypix = []
        for i in range(count):
            x, y = self.render_tools._coordinate_to_pixel(self.gas.positions[i])
            self.xpix.append(x)
            self.ypix.append(y)

    def _deposit_nearest(self, x, y, f, w, weighted, weights):
        weighted[y, x] += f * w
        weights[y, x] += w

    def _deposit_gaussian(self, h, x, y, f, w, weighted, weights):
        """
        Deposit a quantity onto the image using a Gaussian kernel.
        Includes extensive sanity checks for debugging.
        """
        if not np.isfinite(h):
            print(f"Invalid gas size h = {h}")
            return
        
        sigma = self.render_tools._pixel_radius(h)
        
        if not np.isfinite(sigma):
            print(f"Sigma is NaN/Inf. h={h}")
            return
        
        if sigma <= 0:
            print(f"Sigma <=0. sigma = {sigma}, h = {h}")
            return

        if not np.isfinite(f):
            print(f"Field value is invalid: {f}")
            return
        
        if not np.isfinite(w):
            print(f"Weight value is invalid: {w}")
            return

        if not (0 <= x < self.scene.buff_size):
            print(f"x outside image: {x}")
            return

        if not (0 <= y < self.scene.buff_size):
            print(f"y outside image: {y}")
            return

        pixel_radius = int(np.ceil(3 * sigma))

        xmin = max(0, x-pixel_radius)
        xmax = min(self.scene.buff_size, x + pixel_radius + 1)

        ymin = max(0, y-pixel_radius)
        ymax = min(self.scene.buff_size, y + pixel_radius + 1)

        #-----------------------------------------------------
        # Build kernel first
        #-----------------------------------------------------

        kernel = []

        for xpix in range(xmin, xmax):
            for ypix in range(ymin, ymax):

                dx = xpix - x
                dy = ypix - y
                r2 = dx * dx + dy * dy

                exponent = -r2 / (2 * sigma * sigma)

                if not np.isfinite(exponent):
                    print(f"Bad exponent = {exponent}")
                    print(f"r2={r2}, sigma={sigma}")
                    continue

                weight = np.exp(exponent)

                if not np.isfinite(weight):
                    print(f"Gaussian weight became Nan/Inf")
                    print(f"exponent={exponent}")
                    continue

                if weight < 1e-4:
                    continue

                kernel.append((xpix, ypix, weight))
        if len(kernel) == 0:
            print(f"No kernel pixels for sigma = {sigma:.2f}")
            return

        #-----------------------------------------------------
        # Normalize kernel
        #-----------------------------------------------------

        kernel_sum = sum(weight for _, _, weight in kernel)

        if not np.isfinite(kernel_sum):
            print("Kernel sum is NaN/Inf")
            return

        if kernel_sum <= 0:
            print("Kernel sum <=0")
            return
        

        for xpix, ypix, weight in kernel:

            weight /= kernel_sum
            contribution = f * w * weight

            if not np.isfinite(contribution):
                print("Contribution became NaN")
                print(f"f={f}")
                print(f"w={w}")
                print(f"weight={weight}")
                continue
        
            weighted[ypix, xpix] += f * w * weight
            weights[ypix, xpix] += w * weight

        # ------------------------------------------------------------
        # Final sanity check
        # ------------------------------------------------------------
        if np.isnan(weighted).any():
            idx = np.argwhere(np.isnan(weighted))[0]
            print(f"NaN introduced into weighted at pixel {tuple(idx)}")

        if np.isnan(weights).any():
            idx = np.argwhere(np.isnan(weights))[0]
            print(f"NaN introduced into weights at pixel {tuple(idx)}")

    def _deposit_sph(self, r, x, y, f, w, weighted, weights):

        h = self.render_tools._pixel_radius(r)

        radius = int(np.ceil(2 * h))

        xmin = max(0, x-radius)
        xmax = min(self.scene.buff_size, x+radius+1)

        ymin = max(0, y-radius)
        ymax = min(self.scene.buff_size, y+radius+1)

        for xpix in range(xmin, xmax):
            for ypix in range(ymin, ymax):

                dx = xpix - x
                dy = ypix - y

                dist = np.sqrt(dx*dx + dy*dy)

                q = dist / h

                if q >= 2:
                    continue

                if q < 1:
                    kernel = 1 - 1.5*q*q + 0.75*q*q*q
                else:
                    kernel = 0.25*(2-q)**3

                weighted[ypix, xpix] += f * w * kernel
                weights[ypix, xpix] += w * kernel

    def deposit_field(self, field, weight_field="masses", kernel="nearest"):


        if isinstance(field, str):
            field = getattr(self.gas, field)

        if isinstance(weight_field, str):
            weight_field = getattr(self.gas, weight_field)

        weighted = np.zeros((self.scene.buff_size, self.scene.buff_size), dtype=float)
        weights = np.zeros_like(weighted)
        
        for h, x, y, f, w in zip(self.gas.sizes, self.xpix, self.ypix, field, weight_field):
            if kernel == "nearest":
                self._deposit_nearest(x, y, f, w, weighted, weights)
            elif kernel == "gaussian":
                self._deposit_gaussian(h, x, y, f, w, weighted, weights)
            elif kernel == "sph":
                self._deposit_sph(h, x, y, f, w, weighted, weights)

        grid = np.divide(weighted, weights, out=np.zeros_like(weighted), where=weights>0)
        
        return grid, weighted, weights

    def set_background(cmap):
        """
        Return a copy of the colormap with masked values rendered
        using its minimum color.
        """

        if isinstance(cmap, str):
            cmap = plt.get_cmap(cmap)

        cmap = cmap.copy()

        background = cmap(0.0)      # first color of the colormap
        cmap.set_bad(background)

        return cmap

    def radial_velocity(self):

        center = self.scene.center 

        r = self.gas.positions - center
        rmag = np.linalg.norm(r, axis=1)

        bulk_velocity = np.average(self.gas.velocities, axis=0, weights=self.gas.masses)

        v = self.gas.velocities - bulk_velocity

        vr = np.zeros_like(rmag)
        mask = rmag > 0
        vr[mask] = np.sum(r[mask] * v[mask], axis=1) / rmag[mask]

        return vr        


        

"""
        self.snap = snap
        sim.binary_reader.read_sink_snap(self.snap)
        self.snapshot_time = find_cosmic_time(sim.binary_reader.time)
        
        self.galaxy_ind = galaxy_ind
        snapshot = self.sim.snapshot(self.snap)
        self.Galaxies = snapshot.galaxies
        Galaxy = self.Galaxies.find_subhalo(self.galaxy_ind)
        self.center = Galaxy.position
        self.radius = Galaxy.radius * span
        self.bh_ids = Galaxy.bh_particle_ids

        self.render_tools = RenderTools(self)

        self.yt_utils = snapshot.yt

        region = self.yt_utils.sphere(self.center, self.radius, units="code_length")
        
        self.pos4 = region.data[("PartType4", "Coordinates")].v
        self.mass4 = region.data[("PartType4", "Masses")].in_units("Msun").value
        self.id4 = region.data[("PartType4", "ParticleIDs")].v

        self.id4_lookup = {
            pid: i
            for i, pid in enumerate(self.id4)
        }

        self.pos5 = region.data[("PartType5", "Coordinates")].v
        self.id5 = region.data[("PartType5", "ParticleIDs")].v
        self.mass5 = region.data[("PartType5", "Masses")].in_units("Msun").value

        self.id5_lookup = {
            pid: i
            for i, pid in enumerate(self.id5)
        }
    
    def star_renderer(self):
        return StellarRenderer(self)



    for pid in self.popii_particles:

    if pid in self.scene.id4:
        i = self.scene.id4_lookup[pid]
        p4 = self.scene.pos4[i][0]
        m4 = self.scene.mass4[i]

        xpix, ypix = self.render_tools._coordinate_to_pixel(p4)

        meta = self.popii_particles[pid]["meta"]
        formation_time = find_cosmic_time(meta["FormationTime"])
        age = (self.snapshot_time - formation_time) / 1e6 # in Myr

        color, lum = self.color_by_age(age)

        L = lum * m4

        R[ypix, xpix] += L * color[0]
        G[ypix, xpix] += L * color[1]
        B[ypix, xpix] += L * color[2]
    
    if pid in self.scene.id5:
        i = self.scene.id5_lookup[pid]
        p5 = self.scene.pos5[i][0]
        m5 = self.scene.mass5[i]

        xpix, ypix = self.render_tools._coordinate_to_pixel(p5)

        meta = self.popii_particles[pid]["meta"]
        formation_time = find_cosmic_time(meta["FormationTime"])
        age = (self.scene.snapshot_time - formation_time) / 1e6 # in Myr

        color, lum = self.color_by_age(age)

        L = lum * m5

        R[ypix, xpix] += L * color[0]
        G[ypix, xpix] += L * color[1]
        B[ypix, xpix] += L * color[2]

RGB = self.render_tools._post_processing(R, G, B)




    def render(self):

        number_density = np.zeros(self.scene.buff_size)
        


        rho = np.load(f"{self.sim.label}/frb_density_156_{self.galaxy_ind}.npy")
        print(np.shape(rho))
        T = np.load(f"{self.sim.label}/frb_temp_156_{self.galaxy_ind}.npy")
        Z = np.load(f"{self.sim.label}/frb_metallicity_156_{self.galaxy_ind}.npy")

        I = np.arcsinh(rho / np.percentile(rho, 99))
        I /= I.max()

        logT = np.log10(T)
        Tmin = np.min(logT)
        Tmax = np.max(logT)
        T = (logT - Tmin) / (Tmax - Tmin)
        T = np.clip(T, 0, 1)

        cmap = plt.get_cmap("coolwarm")

        color = cmap(T)
        
        color = color[:, :, :3]
        
        RGB = color * I[:, :, None]
        print(np.shape(RGB))
        #S = normalize(metallicity)
        #RGB = RGB*(0.5+0.5*S[:,:,None])

        return RGB
"""