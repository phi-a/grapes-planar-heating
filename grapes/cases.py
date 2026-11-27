"""Reference entry cases."""
from __future__ import annotations

import numpy as np

from dataclasses import replace

from .atmosphere import mars_pathfinder_1997
from .bodies import EARTH, MARS, TITAN
from .frames import inertial_to_relative
from .vehicle import Vehicle

# --- Stardust sample return capsule, Earth entry 2006-01-15 -----------------
# Sources:
#   [D00] Desai, Mitcheltree, Cheatwood, "Sample Returns Missions in the
#         Coming Decade", IAF-00-Q.2.04 (NASA NTRS 20000109794)
#   [D08] Desai, Qualls, "Stardust Entry Reconstruction", AIAA-2008-1198
#         (NASA NTRS 20080008567)
#   [M97] Mitcheltree, Wilmoth, Cheatwood, Brauckmann, Greene, "Aerodynamics
#         of Stardust Sample Return Capsule", AIAA 97-2304 (NTRS 20040105538)
STARDUST = Vehicle(
    name="Stardust SRC",
    mass=46.0,  # kg [D00]
    ref_area=0.51887,  # m^2, frontal area, 0.8128 m diameter [M97]
    cd=1.49,  # zero-angle axial coefficient, 1.482-1.499 from Mach 40 to 12 [M97 Table 4]
    nose_radius=0.2286,  # m [M97]
)

# Space-fixed state at entry interface (radius 6503.14 km) [D08].
STARDUST_INERTIAL = dict(
    h0=135e3,  # m; that radius is ~135 km above the surface at this latitude [D08 Fig. 8]
    v=12.9e3,  # m/s
    gamma=np.radians(-8.2),
    azimuth=np.radians(102.9),  # entry azimuth, east-southeast [D08]
    latitude=np.radians(42.0),  # estimated: landing latitude 40.37 N plus the ground track
)

# Flight values to check against [D08].
STARDUST_FLIGHT = dict(
    peak_g=32.89,  # best-estimated trajectory; pre-entry prediction 32.86
    drogue_time=137.9,  # s after entry interface
    drogue_altitude=31.03e3,  # m MSL
    drogue_mach=1.23,
    peak_q_design=1200.0,  # W/cm^2, convective + radiative, pre-flight [D00]
)


def stardust_convert(v, gamma):
    """Space-fixed (v, gamma) at Stardust's entry point -> atmosphere-relative."""
    s = STARDUST_INERTIAL
    return inertial_to_relative(EARTH, s["h0"], v, gamma, s["latitude"], s["azimuth"])[:2]


def _entry(body, vehicle, inertial):
    """Arguments for simulate(): atmosphere-relative state plus rotation terms."""
    v, gamma, azimuth = inertial_to_relative(
        body, inertial["h0"], inertial["v"], inertial["gamma"],
        inertial["latitude"], inertial["azimuth"])
    return dict(body=body, vehicle=vehicle, h0=inertial["h0"], v0=float(v),
                gamma0=float(gamma), latitude=inertial["latitude"], azimuth=float(azimuth))


# Ready to pass to simulate().
STARDUST_ENTRY = _entry(EARTH, STARDUST, STARDUST_INERTIAL)


# --- Mars Pathfinder, Mars entry 1997-07-04 ---------------------------------
# Sources:
#   [B97] Braun, Spencer, Kallemeyn, Vaughan, "Mars Pathfinder Atmospheric
#         Entry Navigation Operations", AIAA 97-3663 (NASA TM-112885)
#   [G95] Gnoffo, Weilmuenster, Braun, Cruz, "Effects of Sonic Line Transition
#         on Aerothermodynamics of the Mars Pathfinder Probe", AIAA 1995
PATHFINDER = Vehicle(
    name="Mars Pathfinder",
    mass=585.3,  # kg [B97]
    ref_area=np.pi * (2.65 / 2) ** 2,  # m^2, 2.65 m diameter [G95 Fig. 1]
    cd=585.3 / (62.4 * np.pi * (2.65 / 2) ** 2),  # 1.70, from ballistic coefficient 62.4 kg/m^2 [B97]
    nose_radius=0.6625,  # m [G95 Fig. 1]
)

_PATHFINDER_SITE_RADIUS = 3389.72e3  # m, landing site; the atmosphere table is measured from here

# Mars as Pathfinder met it: measured atmosphere, altitudes above the landing site.
MARS_AT_PATHFINDER = replace(MARS, radius=_PATHFINDER_SITE_RADIUS,
                             atmosphere=mars_pathfinder_1997())

# Best-estimated space-fixed state at entry interface [B97 Table 1].
PATHFINDER_INERTIAL = dict(
    h0=3522.2e3 - _PATHFINDER_SITE_RADIUS,  # 132.5 km
    v=7.264e3,
    gamma=np.radians(-14.060),
    azimuth=np.radians(253.148),
    latitude=np.radians(22.630),
)

# The same table's atmosphere-relative state: an independent check on frames.py.
PATHFINDER_RELATIVE_PUBLISHED = dict(v=7.479e3, gamma_deg=-13.649, azimuth_deg=253.675)

# Flight values to check against [B97].
PATHFINDER_FLIGHT = dict(
    peak_g=16.0,  # "16 g", about 78 s after entry
    peak_g_time=78.0,
    peak_q=100.0,  # W/cm^2, "approximately 100", about 70 s after entry
    peak_q_time=70.0,
    chute_time=171.0,  # s after entry
    chute_altitude=9e3,  # m, approximate [B97 Fig. 7]
    chute_speed=375.0,  # m/s [B97 Fig. 7]
)


def pathfinder_convert(v, gamma):
    """Space-fixed (v, gamma) at Pathfinder's entry point -> atmosphere-relative."""
    s = PATHFINDER_INERTIAL
    return inertial_to_relative(MARS_AT_PATHFINDER, s["h0"], v, gamma,
                                s["latitude"], s["azimuth"])[:2]


PATHFINDER_ENTRY = _entry(MARS_AT_PATHFINDER, PATHFINDER, PATHFINDER_INERTIAL)


# --- Huygens, Titan entry 2005-01-14 -----------------------------------------
# Sources:
#   [S07] Striepe, Blanchard, Kirsch, Fowler, "Huygens Titan Probe Trajectory
#         Reconstruction Using Traditional Methods and the Program to Optimize
#         Simulated Trajectories II", AAS 07-226 (NASA NTRS 20070010747)
#   [J06] Johnston, Hollis, Sutton, "Radiative Heating Methodology for the
#         Huygens Probe", AIAA 2006-3426 (NASA NTRS 20080008445)
#   [K05] Kazeminejad, "Methodology Development for the Reconstruction of the
#         ESA Huygens Probe Entry and Descent Trajectory", PhD thesis,
#         University of Graz, 2005
HUYGENS = Vehicle(
    name="Huygens",
    mass=320.0,  # kg, entry configuration [K05 Table 3.1]; held constant (ablation loses 9.7 kg)
    ref_area=5.73,  # m^2, heat shield reference area [S07, K05 Table 3.1]
    cd=1.54,  # hypersonic axial force coefficient, about 1.53-1.56 [S07 Fig. 16]
    nose_radius=1.25,  # m [J06]
)

# JPL's last estimate of the space-fixed state at entry interface [S07 Table 1].
_r = np.array([-3785.053, 366.623, -568.429]) * 1e3  # m
_vel = np.array([5.704492, 1.918924, 0.390306]) * 1e3  # m/s
HUYGENS_INERTIAL = dict(
    h0=float(np.linalg.norm(_r) - TITAN.radius),  # 1270 km
    v=float(np.linalg.norm(_vel)),  # 6031 m/s
    gamma=float(np.arcsin(_r @ _vel / (np.linalg.norm(_r) * np.linalg.norm(_vel)))),  # -65.6 deg
)

# Titan turns so slowly (about 17 m/s at entry altitude) that rotation is
# left out: the space-fixed state is used directly.
HUYGENS_ENTRY = dict(
    body=TITAN,
    vehicle=HUYGENS,
    h0=HUYGENS_INERTIAL["h0"],
    v0=HUYGENS_INERTIAL["v"],
    gamma0=HUYGENS_INERTIAL["gamma"],
)

# Flight values. Times are seconds after entry interface (09:05:52.523 UTC [S07]).
HUYGENS_FLIGHT = dict(
    peak_decel=121.0,  # m/s^2, read from [S07 Fig. 14]
    peak_decel_time=198.5,  # s, read from [S07 Fig. 14]
    late_time=267.7,  # s: last record of the HASI entry profile (09:10:20.188 UTC)
    late_altitude=157.318e3,  # m, altitude of that record
    pilot_chute_time=273.0,  # s, approximate [S07 Fig. 19]
)
