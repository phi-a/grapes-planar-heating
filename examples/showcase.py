"""Portfolio figures from the 2-D model. Run: python examples/showcase.py

  stardust_entry.png      - Stardust at Earth: trajectory, heating, deceleration
  stardust_corridor.png   - why Stardust came in at -8.2 deg: the entry corridor
  pathfinder_entry.png    - Pathfinder at Mars, same three panels
  pathfinder_corridor.png - Pathfinder's place on the Mars heating map
  huygens_entry.png       - Huygens at Titan, same three panels
  three_worlds.png        - convective vs radiative peak heating, all three
  size_effect.png         - radiative share of peak heating vs nose radius
"""
from dataclasses import replace
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap

from grapes import TITAN, simulate
from grapes.cases import (HUYGENS_ENTRY, HUYGENS_FLIGHT, PATHFINDER_ENTRY, PATHFINDER_FLIGHT,
                          PATHFINDER_INERTIAL, STARDUST_ENTRY, STARDUST_FLIGHT,
                          STARDUST_INERTIAL, pathfinder_convert, stardust_convert)
from grapes.heating import HUYGENS_DESIGN_POINTS, TITAN_HUYGENS_BOLTZMANN
from grapes.sweep import entry_corridor, skip_out_boundary

OUT = Path(__file__).resolve().parents[1] / "docs" / "figures"

# Palette: one blue ramp for magnitude; blue/orange for convective/radiative.
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
BLUE_RAMP = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "font.family": "DejaVu Sans", "font.size": 10,
    "text.color": INK, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.edgecolor": GRID, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "axes.titlesize": 11, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "lines.linewidth": 2.0, "lines.solid_capstyle": "round",
})


def _rotation(entry):
    return dict(latitude=entry["latitude"], azimuth=entry["azimuth"])


def entry_figure(entry, filename, title, t_end, late_point, g_ref, q_ref=None, g_limit=None,
                 q_max=None, g_max=None, h_range=(0, 140), t_start=0.0, rad_upper=None,
                 design=None, late_offset=(14, 4), late_compact=False):
    """Three panels: trajectory, heating, deceleration.

    late_point: (time s, flight altitude m, label) for a late-entry check.
    q_ref / g_ref: (value, label) reference line / flight value.
    rad_upper: (radiative model, label) for an upper radiative estimate.
    design: published design-study points to mark on the heating panel.
    """
    tr = simulate(**entry, t_end=t_end, dt_out=0.1)
    g = tr.decel / 9.80665

    fig, ax = plt.subplots(1, 3, figsize=(13, 4.2), constrained_layout=True)
    fig.suptitle(title, x=0.01, ha="left", fontsize=14, fontweight="bold")

    a = ax[0]
    a.plot(tr.v / 1e3, tr.h / 1e3, color=BLUE)
    i = tr.q_total.argmax()
    a.plot(tr.v[i] / 1e3, tr.h[i] / 1e3, "o", color=BLUE, ms=8, mec=SURFACE, mew=2)
    a.annotate(f"peak heating\n{tr.h[i] / 1e3:.0f} km", (tr.v[i] / 1e3, tr.h[i] / 1e3),
               xytext=(-10, 10), textcoords="offset points", ha="right", color=INK2)
    t_late, h_late, what = late_point
    k = np.argmin(abs(tr.t - t_late))
    a.plot(tr.v[k] / 1e3, h_late / 1e3, "o", color=SURFACE, ms=8, mec=INK, mew=2)
    sep = " " if late_compact else "\n"
    a.annotate(f"at {tr.t[k]:.0f} s: grapes {tr.h[k] / 1e3:.1f} km,{sep}"
               f"flight {h_late / 1e3:.1f} km" + (f" ({what})" if what else ""),
               (tr.v[k] / 1e3, h_late / 1e3), xytext=late_offset, textcoords="offset points",
               va="center", color=INK)
    a.set(title="Trajectory", xlabel="Speed through the air (km/s)", ylabel="Altitude (km)",
          ylim=h_range)

    a = ax[1]
    # Where the radiative model has no data, leave a gap rather than draw a total.
    total = np.where(np.isnan(tr.q_rad), np.nan, tr.q_total)
    a.plot(tr.t, total / 1e4, color=INK, label=f"total, peak {np.nanmax(total) / 1e4:.0f}")
    a.plot(tr.t, tr.q_conv / 1e4, color=BLUE,
           label=f"convective, peak {tr.peak_q_conv / 1e4:.0f}")
    a.plot(tr.t, tr.q_rad / 1e4, color=ORANGE,
           label=f"radiative, peak {tr.peak_q_rad / 1e4:.0f}")
    q_top = q_max or 1.1 * max(q_ref[0] if q_ref else 0, tr.peak_q_total / 1e4)
    if rad_upper:
        model, label = rad_upper
        up = simulate(**{**entry, "body": replace(entry["body"], radiative=model)},
                      t_end=t_end, dt_out=0.1)
        a.plot(up.t, up.q_rad / 1e4, color=ORANGE, ls=(0, (4, 2)),
               label=f"{label}, peak {up.peak_q_rad / 1e4:.0f}")
    if design:
        a.plot(design["t"], design["q_conv"], "o", color=SURFACE, mec=BLUE, mew=2, ms=7,
               label="design, convective")
        a.plot(design["t"], design["q_rad_nonboltzmann"], "o", color=SURFACE, mec=ORANGE,
               mew=2, ms=7, label="design, radiative")
        a.plot(design["t"], design["q_rad_boltzmann"], "s", color=SURFACE, mec=ORANGE,
               mew=2, ms=7, label="design, radiative (older)")
    if q_ref:
        q_val, q_label = q_ref
        a.axhline(q_val, color=INK2, lw=1, ls=(0, (4, 3)))
        a.text(t_start, q_val + 0.015 * q_top, q_label, ha="left", va="bottom", color=INK2)
    a.set(title="Heating at the nose (W/cm²)", xlabel="Time from entry (s)", ylim=(0, q_top),
          xlim=(t_start - 0.03 * (t_end - t_start), t_end + 0.03 * (t_end - t_start)))
    if design:
        a.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.56, 1.02), fontsize=8,
                 handlelength=1.6)
    else:
        a.legend(frameon=False, loc="upper right", bbox_to_anchor=(1.0, 0.66))

    a = ax[2]
    a.plot(tr.t, g, color=BLUE)
    g_val, g_label = g_ref
    g_top = g_max or 1.25 * max(g_val, g.max(), g_limit or 0)
    if g_limit:
        a.axhline(g_limit, color=INK2, lw=1, ls=(0, (4, 3)))
        a.text(t_end - 2, g_limit - 0.02 * g_top, f"design limit {g_limit:.0f} g",
               ha="right", va="top", color=INK2)
    j = g.argmax()
    a.plot(tr.t[j], g[j], "o", color=BLUE, ms=8, mec=SURFACE, mew=2)
    a.annotate(f"grapes {g[j]:.1f} g\n{g_label}", (tr.t[j], g[j]),
               xytext=(12, -4), textcoords="offset points", va="top", color=INK)
    a.set(title="Deceleration (Earth g)", xlabel="Time from entry (s)", ylim=(0, g_top),
          xlim=(t_start - 0.03 * (t_end - t_start), t_end + 0.03 * (t_end - t_start)))

    fig.savefig(OUT / filename, dpi=180)
    plt.close(fig)
    return tr


def corridor_figure(entry, inertial, convert, filename, title, speeds, angles_deg, levels,
                    point_label, g_limit=None, q_limit=None, g_label_at=None,
                    q_label_at=None):
    """Peak total heating over space-fixed entry speed and angle."""
    rot = _rotation(entry)
    m = entry_corridor(entry["body"], entry["vehicle"], entry["h0"], speeds,
                       np.radians(angles_deg), convert=convert, **rot)
    edge = np.degrees(skip_out_boundary(
        entry["body"], entry["vehicle"], entry["h0"], speeds,
        np.radians(angles_deg[0]), np.radians(angles_deg[-1]), convert=convert, **rot))
    q, gload, V = m.peak_q_total / 1e4, m.peak_g, speeds / 1e3
    top = angles_deg[0]

    fig, a = plt.subplots(figsize=(9, 6), constrained_layout=True)
    cmap = ListedColormap(BLUE_RAMP)
    cs = a.contourf(V, angles_deg, q, levels=levels, cmap=cmap,
                    norm=BoundaryNorm(levels, cmap.N), extend="max")
    a.contour(V, angles_deg, q, levels=levels[1:], colors=SURFACE, linewidths=1.2)
    cb = fig.colorbar(cs, ax=a, pad=0.02, ticks=levels)
    cb.set_label("Peak heating at the nose, convective + radiative (W/cm²)", color=INK2)
    cb.outline.set_visible(False)

    if g_limit:
        cg = a.contour(V, angles_deg, gload, levels=[g_limit], colors=[INK], linewidths=2)
        a.clabel(cg, fmt=f"{g_limit:.0f} g limit", inline=True, fontsize=10,
                 manual=[g_label_at])
    if q_limit:
        a.contour(V, angles_deg, q, levels=[q_limit], colors=[ORANGE], linewidths=2.5)
        a.text(*q_label_at, f"{q_limit:.0f} W/cm²\nheating limit", color=SURFACE,
               ha="left", va="center")

    if not np.all(np.isnan(edge)):
        a.fill_between(V, np.where(np.isnan(edge), top, edge), top, color=GRID, zorder=3,
                       lw=0)
        a.plot(V, edge, color=SURFACE, lw=2, zorder=4)
        a.text(V[-1] - 0.02 * (V[-1] - V[0]), top - 0.035 * (top - angles_deg[-1]),
               "skips back out to space", ha="right", va="center", color=INK2, zorder=5)

    v_pt, g_pt = inertial["v"] / 1e3, np.degrees(inertial["gamma"])
    a.plot(v_pt, g_pt, "o", color=SURFACE, ms=11, mec=INK, mew=2, zorder=6)
    a.annotate(point_label, (v_pt, g_pt), xytext=(-12, 2), textcoords="offset points",
               ha="right", va="center", color=INK, fontweight="bold", zorder=6)

    a.set(xlabel="Entry speed, space-fixed (km/s)",
          ylabel="Entry flight-path angle, space-fixed (deg)",
          xlim=(V[0], V[-1]), ylim=(angles_deg[-1], top))
    a.set_title(title, fontsize=12)
    a.grid(False)
    fig.savefig(OUT / filename, dpi=180)
    plt.close(fig)
    return m


def three_worlds(runs, filename):
    """Peak convective and radiative heating side by side for each mission.

    runs: list of (label, trajectory, optional upper-radiative trajectory).
    """
    fig, ax = plt.subplots(1, len(runs), figsize=(11, 3.9), constrained_layout=True)
    fig.suptitle("Where the heat comes from: convective vs radiative at three worlds",
                 x=0.01, ha="left", fontsize=14, fontweight="bold")
    for a, (label, tr, upper) in zip(ax, runs):
        conv, rad = tr.peak_q_conv / 1e4, tr.peak_q_rad / 1e4
        a.set_axisbelow(True)
        bars = a.bar([0, 1], [conv, rad], width=0.62, color=[BLUE, ORANGE])
        top = max(conv, rad)
        shares = [rad / (conv + rad)]
        if upper is not None:
            hi = upper.peak_q_rad / 1e4
            a.plot([1, 1], [rad, hi], color=INK2, lw=1.5)
            a.plot([0.86, 1.14], [hi, hi], color=INK2, lw=1.5)
            a.text(1.2, hi, f"{hi:.0f} with the\nolder model", va="center", color=INK2,
                   fontsize=9)
            top = max(top, hi)
            shares.append(hi / (conv + hi))
        for b, val, dx in zip(bars, (conv, rad), (0.0, -0.17 if upper is not None else 0.0)):
            a.text(b.get_x() + b.get_width() / 2 + dx, val + 0.02 * top, f"{val:.0f}",
                   ha="center", va="bottom", color=INK)
        share = "–".join(f"{100 * x:.0f}" for x in shares)
        a.set_title(f"{label}\nradiative share: {share}%", fontsize=10.5)
        a.set_xticks([0, 1], ["convective", "radiative"])
        a.set_xlim(-0.6, 2.2 if upper is not None else 1.6)
        a.set_ylim(0, 1.2 * top)
        a.grid(axis="x", visible=False)
        a.tick_params(axis="x", length=0)
    ax[0].set_ylabel("Peak heating at the nose (W/cm²)")
    fig.savefig(OUT / filename, dpi=180)
    plt.close(fig)


def radiative_share(entry, nose_radius, t_end):
    """Peak radiative / (peak radiative + peak convective), trajectory unchanged."""
    vehicle = replace(entry["vehicle"], nose_radius=nose_radius)
    tr = simulate(**{**entry, "vehicle": vehicle}, t_end=t_end, dt_out=0.1)
    return tr.peak_q_rad / (tr.peak_q_rad + tr.peak_q_conv)


def size_figure(cases, titan, filename):
    """Radiative share of peak heating against nose radius.

    cases: (label, color, entry, t_end), one line each, drawn over the nose
    radii the body's radiative fit covers and dashed down to the flown radius.
    titan: (label, radius, low share, high share), a single flown point.
    """
    fig, a = plt.subplots(figsize=(8.2, 4.6), constrained_layout=True)
    a.axhline(50, color=INK2, lw=1, ls=(0, (4, 3)))
    a.text(85, 51.5, "radiation exceeds convection", color=INK2, va="bottom", ha="right")
    for label, color, entry, t_end in cases:
        flown = entry["vehicle"].nose_radius
        lo, hi = entry["body"].radiative.valid["nose_radius"]
        radii = np.geomspace(lo, hi, 40)
        share = np.array([100 * radiative_share(entry, r, t_end) for r in radii])
        a.plot(radii, share, color=color)
        if flown < lo:
            a.plot([flown, lo], [100 * radiative_share(entry, flown, t_end), share[0]],
                   color=color, ls=(0, (2, 2)))
        f = 100 * radiative_share(entry, flown, t_end)
        a.plot(flown, f, "o", color=color, ms=8, mec=SURFACE, mew=2)
        cross = np.interp(50, share, radii)
        a.plot(cross, 50, "o", color=SURFACE, ms=7, mec=color, mew=2)
        a.annotate(f"{cross:.1f} m", (cross, 50), xytext=(8, -7), textcoords="offset points",
                   ha="left", va="top", color=INK)
        a.annotate(label, (radii[-1], share[-1]), xytext=(6, 0), textcoords="offset points",
                   va="center", color=INK)
        print(f"{label.splitlines()[0]}: {f:.0f}% as flown ({flown:.2f} m), 50% at {cross:.1f} m")
    label, radius, low, high = titan
    a.plot([radius, radius], [low, high], color=AQUA, lw=2)
    a.plot([radius, radius], [low, high], "o", color=AQUA, ms=8, mec=SURFACE, mew=2)
    a.annotate(label, (radius, high), xytext=(-9, 0), textcoords="offset points", va="center",
               ha="right", color=INK)
    a.set_xscale("log")
    a.set_xticks([0.1, 0.3, 1, 3, 10, 30], ["0.1", "0.3", "1", "3", "10", "30"])
    a.set(xlim=(0.1, 90), ylim=(0, 80), xlabel="Nose radius (m)",
          ylabel="Radiative share of peak heating (%)")
    a.set_title("Radiative share of peak heating grows with nose radius", fontsize=12)
    fig.savefig(OUT / filename, dpi=180)
    plt.close(fig)


def report(name, tr):
    print(f"{name}: peak {tr.peak_q_conv / 1e4:.0f} convective + {tr.peak_q_rad / 1e4:.0f} "
          f"radiative, {tr.peak_q_total / 1e4:.0f} total W/cm2; peak {tr.peak_decel_g:.1f} g")


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    SF, PF = STARDUST_FLIGHT, PATHFINDER_FLIGHT

    sd = entry_figure(
        STARDUST_ENTRY, "stardust_entry.png",
        "Stardust capsule entering Earth's atmosphere at 12.9 km/s", t_end=140,
        late_point=(SF["drogue_time"], SF["drogue_altitude"], "drogue out"),
        q_ref=(SF["peak_q_design"], "pre-flight design peak 1200"),
        g_ref=(SF["peak_g"], f"flight {SF['peak_g']:.1f} g"), g_limit=40.0,
        q_max=1300, g_max=44)
    report("Stardust", sd)

    pf = entry_figure(
        PATHFINDER_ENTRY, "pathfinder_entry.png",
        "Mars Pathfinder entering Mars's atmosphere at 7.3 km/s", t_end=PF["chute_time"],
        late_point=(PF["chute_time"], PF["chute_altitude"], "parachute out"),
        q_ref=(PF["peak_q"], "reported ≈ 100"),
        g_ref=(PF["peak_g"], f"flight {PF['peak_g']:.0f} g"),
        q_max=140, g_max=20)
    report("Pathfinder", pf)

    HF = HUYGENS_FLIGHT
    hy = entry_figure(
        HUYGENS_ENTRY, "huygens_entry.png",
        "Huygens entering Titan's atmosphere at 6.0 km/s", t_end=HF["late_time"],
        late_point=(HF["late_time"], HF["late_altitude"], ""),
        g_ref=(HF["peak_decel"] / 9.80665, f"flight {HF['peak_decel'] / 9.80665:.1f} g"),
        q_max=95, g_max=15, h_range=(0, 1300), t_start=140.0, late_offset=(14, -16),
        late_compact=True,
        rad_upper=(TITAN_HUYGENS_BOLTZMANN, "radiative (older)"),
        design=HUYGENS_DESIGN_POINTS)
    report("Huygens", hy)
    hy_upper = simulate(**{**HUYGENS_ENTRY, "body": replace(TITAN, radiative=TITAN_HUYGENS_BOLTZMANN)},
                        t_end=HF["late_time"], dt_out=0.1)
    three_worlds([("Earth: Stardust, 12.9 km/s", sd, None),
                  ("Mars: Pathfinder, 7.3 km/s", pf, None),
                  ("Titan: Huygens, 6.0 km/s", hy, hy_upper)], "three_worlds.png")
    size_figure(
        [("Earth\nStardust's entry", BLUE, STARDUST_ENTRY, 140),
         ("Mars\nPathfinder's entry", ORANGE, PATHFINDER_ENTRY, PF["chute_time"])],
        ("Titan, Huygens\n(one size only)", HUYGENS_ENTRY["vehicle"].nose_radius,
         100 * hy.peak_q_rad / (hy.peak_q_rad + hy.peak_q_conv),
         100 * hy_upper.peak_q_rad / (hy_upper.peak_q_rad + hy.peak_q_conv)),
        "size_effect.png")

    corridor_figure(
        STARDUST_ENTRY, STARDUST_INERTIAL, stardust_convert, "stardust_corridor.png",
        "Stardust's entry corridor: too shallow skips out, too steep exceeds 40 g",
        speeds=np.linspace(8e3, 16e3, 41), angles_deg=np.linspace(-4.0, -16.0, 49),
        levels=[0, 200, 400, 600, 800, 1200, 1600, 2400],
        point_label="Stardust\n12.9 km/s, −8.2°", g_limit=40.0, q_limit=1200.0,
        g_label_at=(10.0, -11.5), q_label_at=(12.62, -14.6))

    corridor_figure(
        PATHFINDER_ENTRY, PATHFINDER_INERTIAL, pathfinder_convert, "pathfinder_corridor.png",
        "Mars entry heating for the Pathfinder capsule",
        speeds=np.linspace(4e3, 9e3, 41), angles_deg=np.linspace(-8.0, -24.0, 49),
        levels=[0, 25, 50, 75, 100, 150, 200, 300],
        point_label="Pathfinder\n7.26 km/s, −14.1°")
    print("Figures in", OUT)
