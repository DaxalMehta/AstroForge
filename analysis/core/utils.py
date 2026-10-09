from .constants import CS
import numpy as np
import os
from scipy import integrate
from scipy.spatial import cKDTree
from astropy.cosmology import Planck15 as cosmo

def file_range(base, mode=1):
    files = os.listdir(base)
    if mode == 1:
        snapshot_files = [file for file in files if file.startswith("sink_snap_")]
    elif mode == 2:
        snapshot_files = [file for file in files if file.startswith("groups_")]
    elif mode == 3:
        snapshot_files = [file for file in files if file.startswith("Galaxies")]
    else:
        raise ValueError("Incorrect mode entered.")
    if not snapshot_files:
        return None, None
    if mode == 1 or mode == 2:
        snapshot_numbers = [int(file.split("_")[-1]) for file in snapshot_files]
    elif mode == 3:
        snapshot_numbers = [int(file.split(".")[0].split("_")[-1]) for file in snapshot_files]
    start_snap = min(snapshot_numbers)
    end_snap = max(snapshot_numbers) + 1
    return start_snap, end_snap
    
def Type(t):
    if t == 0:
        return "PopIII"
    elif t == 2:
        return "PopII"
    elif t == 3:
        return "MBH"

def SNeType(mass):
    mass *= 1e10 / CS.hubble_parameter

    if (mass < 11): return "No SNe"
    if (mass > 11 and mass <= 40): return "TypeII SNe"
    if (mass > 40 and mass <= 140): return "DCBH"
    if (mass > 140 and mass <= 260): return "PISN"
    if (mass > 260): return "DCBH"

def sort_sink_particles(ids, masses, positions, types):

    mask = types == 0
    popiii_ids = ids[mask]
    popiii_mass = masses[mask]
    popiii_pos = positions[mask]
    
    mask = types == 2
    popii_ids = ids[mask]
    popii_mass = masses[mask]
    popii_pos = positions[mask]
    
    mask = types == 3
    bh_ids = ids[mask]
    bh_mass = masses[mask]
    bh_pos = positions[mask]
    
    return popiii_ids, popiii_mass, popiii_pos, popii_ids, popii_mass, popii_pos, bh_ids, bh_mass, bh_pos

def sink_particles_in_galaxy(Pos, Rad, ids, masses, positions, types):
    
    tree = cKDTree(positions)

    inds = tree.query_ball_point(Pos, r=Rad)

    ids, masses, positions, types = ids[inds], masses[inds], positions[inds], types[inds]

    return sort_sink_particles(ids, masses, positions, types)

def histogram_bins(data, num, type="log"):
    if type == 'log':
        return np.logspace(np.log10(np.min(data)), np.log10(np.max(data)), num)
    if type == 'lin':
        return np.linspace(np.min(data), np.max(data), num)

def unwrap(x):
    if isinstance(x, (list, tuple)) and len(x) == 1:
        return x[0]
    if hasattr(x, "size") and getattr(x, "size", None) == 1:
        return x.item()
    return x

def DotheyAlreadyExist(base, filename):
    files = os.listdir(base)
    if filename in files: return True 
    else: return False

def EddingtonRate(M, e=0.1):
    MdotEd = 2.2e-8 * (0.1 / e) * M
    return MdotEd

def EddingtonLuminosity(M):
    MdotEd = EddingtonRate(M)
    MdotEd *= CS.Msun_to_g / CS.yr_to_s
    EdLuminosity = 0.1 * MdotEd * pow(CS.c_cgs, 2)
    return EdLuminosity

def EddingtonGrowth(M0, t, e=0.1):
    t_Edd = 4.5e7 * (e / 0.1)
    M = M0 * np.exp(t / t_Edd)
    return M

def SuperEddingtonGrowth(M0, t, f_edd=1.0, e=0.1):
    t_Edd = 4.5e7 * (e / 0.1)
    M = M0 * np.exp(f_edd * t / t_Edd)
    return M

def RadiationPressure(M, FEdd, mode=1, PhotonTrapping=False):
    L = EddingtonLuminosity(M)
    if not PhotonTrapping:
        if mode == 1: # Force feedback
            p = FEdd * L / CS.c_in_cgs
            return p
        if mode == 2: # Energy feedback
            E = FEdd * L
            return E
    else:
        if mode == 1:
            p = L * (1 + np.log(FEdd)) / CS.c_in_cgs
            return p
        if mode == 2:
            E = L * (1 + np.log(FEdd))
            return E

def FindEddingtonFactor(M0, M, t, e=0.1):
    t_Edd = 4.5e7 * (e / 0.1)
    f_edd = (t_Edd / t) * np.log(M / M0)
    return f_edd

def FindIMFUncerntainty():

    M_char = 20
    Mmin = 0.1
    Mmax = 300.0
    

    def xi_top(M):
        return M**(-1.3) * np.exp(-(M_char / M)**1.6)
    def xi_sal(M):
        return M**(-2.35)
    def xi_lar(M):
        return M**(-1.35)
    
    NTH, err1 = integrate.quad(xi_top, 40, 140, epsabs=0, epsrel=1e-5)
    print(f"NTH = {NTH:.2e}")
    MtotTH, err2 = integrate.quad(lambda m: m * xi_top(m), Mmin, Mmax, epsabs=0, epsrel=1e-5)

    etaTH = NTH / MtotTH

    NS1, _ = integrate.quad(xi_sal, 40, 140, epsabs=0, epsrel=1e-5)
    NS2, _ = integrate.quad(xi_sal, 260, Mmax, epsabs=0, epsrel=1e-5)

    NS = NS1 + NS2
    print(f"NS = {NS:.2e}, NS1 = {NS1:.2e}, NS2 = {NS2:.2e}")
    MtotS, _ = integrate.quad(xi_lar, Mmin, Mmax, epsabs=0, epsrel=1e-5)
    etaS = NS / MtotS

    R = etaTH / etaS
    print(f"Top-Heavy to Salpeter IMF ratio: {R:.2f}")
    return R

def BondiHoyleAccretion(M, rho, cs, v): # in Msun, g/cm**3, km/s, km/s
    cs *= 1e5
    v *= 1e5
    BondiRadius = M * CS.Msun_to_g * CS.G_in_cgs / (cs**2 + v**2)
    Mdot = 4 * np.pi * rho * BondiRadius**2 * np.sqrt((1.12 * cs)**2 + v**2)
    Mdot = Mdot * CS.yr_to_s / CS.Msun_to_g
    return Mdot

def FreeFallTimescale(mass, radius, mode=1): # in Msun, pc
    mass *= CS.Msun_to_g
    radius *= CS.pc_in_cgs
    if mode == 1: #Sphere infall
        density = 3 * mass / (4 * np.pi * radius**3)
        tff = np.sqrt(3 * np.pi / (32 * CS.G_in_cgs * density))
        return tff
    if mode == 2: #Point mass
        tff = np.pi / 2 * np.sqrt(radius**3 / (2 * CS.G_in_cgs * mass))
        return tff

def PhotonTrappingLuminosity(FEdd, mode=1):
    if mode == 1:
        eta = 1 / 12
        return eta
    if mode == 2:
        eta = 1 - 4 / 3 * (1 / FEdd)**0.5
        return eta
    if mode == 3:
        eta = -1 + 2 / 3 * (1 / FEdd)**0.5 + 0.5 * (FEdd)**0.5
        return eta
    if mode == 4:
        A = pow(0.9963 - 0.9292 * 0.7, -0.5639)
        B = pow(4.627 - 4.445 * 0.7, -0.5524)
        C = pow(827.3 - 718.1 * 0.7, -0.7060)

        eta = A * (0.985/(1/FEdd + B) + 0.015/(1/FEdd + C))
        return eta
    
def CoolingTimeSoundCrossingTimeRation(number_density, temperature, resolution, mu=0.6, f=0.13):
    ## Equation derived from Dalla Vecchia Schaye 2012

    ratio = 2.8e2 * pow(number_density, -1) * (temperature / 10**7.5) * pow(resolution / 100, -1) * pow(mu / 0.6, -3/2) * (f/0.13)

    return ratio

def find_snaps_around_merger(z_merge, snaps, z_snaps):
    z_snaps = np.array(z_snaps)
    snaps = np.array(snaps)
    
    # before
    mask_before = z_snaps >= z_merge
    snap_before = None
    if np.any(mask_before):
        idx = np.argmin(z_snaps[mask_before])
        snap_before = snaps[mask_before][idx]

    # after
    mask_after = z_snaps < z_merge
    snap_after = None
    if np.any(mask_after):
        idx = np.argmin(z_snaps[mask_after])
        snap_after = snaps[mask_after][idx]

    return snap_before, snap_after

def chirp_mass(M1, M2):
    """
    Gives the chirp mass of a given merger event.

    Parameters
    ----------
    M1, M2 : arrays or floats
        BH masses in Msun
    
    Returns
    -------
    M_chirp : same shape
        chirp mass (Msun)
    """

    M1 = np.array(M1)
    M2 = np.array(M2)

    Mtot = M1 + M2
    
    Mchirp = (M1 * M2)**(3/5) / (Mtot)**(1/5)

    return Mchirp

def lisa_detection_probability(M1, M2, z):
    """
    Approximate LISA detection probability (0 → 1)

    Parameters
    ----------
    M1, M2 : arrays or floats
        BH masses in Msun
    z : array or float
        redshift

    Returns
    -------
    P_det : same shape
        detection weight (dimensionless)
    """    

    M1 = np.array(M1)
    M2 = np.array(M2)
    z = np.array(z)

    Mtot = M1 + M2
    
    Mchirp = (M1 * M2)**(3/5) / (Mtot)**(1/5)
    Mchirp_z = Mchirp * (1 + z)

    q = np.minimum(M1, M2) / np.maximum(M1, M2)
    eta = q / (1 + q)**2 

    DL = cosmo.luminosity_distance(z).value # in Mpc

    ### 1. Mass Window (soft edges)
    M_min, M_max = 1e4, 1e7

    w_mass = np.exp(-(np.log10(Mtot/M_min))**-2) * np.exp(-(np.log10(M_max/Mtot))**-2)

    mask_mass = (Mtot > 1e3) * (Mtot < 1e8)
    w_mass *= mask_mass

    ### 2. Chirp mass sensitivity
    w_chirp = np.exp(-(np.log10(Mchirp_z / 1e5))**2 / 2)

    ### 3. Mass ratio weighting
    w_q = eta / 0.25

    ### 4. Distance Suppression
    w_z = 1.0 / (DL**2 + 1e-10)
    w_z /= np.max(w_z)

    P_det = w_mass * w_chirp * w_q * w_z

    if np.max(P_det) > 0:
        P_det /= np.max(P_det)
    
    return P_det

def find_cosmic_time(time):

    redshift = 1 / time - 1
    
    return cosmo.age(redshift).to("yr").value

def sci_notation(x, precision=1):
    if x == 0:
        return "0"

    exponent = int(np.floor(np.log10(abs(x))))
    mantissa = x / 10**exponent

    return rf"{mantissa:.{precision}f}$\times10^{{{exponent}}}$"