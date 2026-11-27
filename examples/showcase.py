"""Figures for the README. Run: python examples/showcase.py

  three_worlds.png        Figure 1  convective vs radiative peak heating
  size_effect.png         Figure 2  radiative share vs nose radius
  validation.png          Figure 3  model minus flight, with Monte Carlo spread
  stardust_entry.png      Figure 4  Stardust entry profile
  pathfinder_entry.png    Figure 5  Pathfinder entry profile
  huygens_entry.png       Figure 6  Huygens entry profile
  stardust_corridor.png   Figure 7  Stardust entry corridor
  pathfinder_corridor.png Figure 8  Pathfinder heating map

Each figure carries its caption at the bottom and no title.
"""
from dataclasses import replace
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap, to_hex, to_rgb
from matplotlib.lines import Line2D

from grapes import TITAN, simulate
from grapes.cases import (HUYGENS_ENTRY, HUYGENS_FLIGHT, PATHFINDER_ENTRY, PATHFINDER_FLIGHT,
                          PATHFINDER_INERTIAL, STARDUST_ENTRY, STARDUST_FLIGHT,
                          STARDUST_INERTIAL, pathfinder_convert, stardust_convert)
from grapes.dispersion import SIGMA, disperse
from grapes.heating import HUYGENS_DESIGN_POINTS, TITAN_HUYGENS_BOLTZMANN
from grapes.sweep import entry_corridor, skip_out_boundary

OUT = Path(__file__).resolve().parents[1] / "docs" / "figures"
G0 = 9.80665

# Palette: ink colors on a white background.
BG, MANILA, INK, GRAPHITE, RULE = "#FFFFFF", "#E9DFCF", "#2D2926", "#65594D", "#B8A995"
SLATE, IRIS, RED, BRASS = "#415571", "#43364F", "#984D39", "#A0783C"
# Roles: Slate = convective, Red pencil = radiative, Ink = total and model lines,
# dashed Graphite/Ink = limits, Brass = Mars, Iris = Titan and flight comparisons.
# Convective is also dashed or hatched, radiative solid or filled, so the two
# stay distinct without color.
WORLD = {"Earth": SLATE, "Mars": BRASS, "Titan": IRIS}


def _ramp(a, b, n):
    a, b = np.array(to_rgb(a)), np.array(to_rgb(b))
    return [to_hex(a + (b - a) * t) for t in np.linspace(0, 1, n)]


HEAT_RAMP = _ramp(MANILA, IRIS, 7)

plt.rcParams.update({
    "figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG,
    "font.family": "DejaVu Sans", "font.size": 9.5,
    "text.color": INK, "axes.labelcolor": GRAPHITE, "xtick.color": GRAPHITE,
    "ytick.color": GRAPHITE, "axes.edgecolor": RULE, "axes.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": "#EFE9DF", "grid.linewidth": 0.9,
    "xtick.major.size": 3, "ytick.major.size": 3,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "legend.frameon": False, "legend.fontsize": 8.5,
    "lines.linewidth": 1.8, "lines.solid_capstyle": "round", "lines.dash_capstyle": "round",
})

DASH = (0, (4, 2.5))


CAPTION_SIZE = 8.6


def _wrap(caption, number, width_in):
    """Greedy wrap by measured text width. The first line starts after the label."""
    probe = plt.figure(figsize=(width_in, 1))
    renderer = probe.canvas.get_renderer()

    def width(text, bold=False):
        t = probe.text(0, 0, text, fontsize=CAPTION_SIZE, fontweight="bold" if bold else None)
        w = t.get_window_extent(renderer).width / probe.dpi
        t.remove()
        return w

    head = f"Figure {number}."
    avail = 0.976 * width_in
    first = avail - width(head + " ", bold=True)
    lines, line = [], ""
    for word in caption.split():
        trial = f"{line} {word}".strip()
        if width(trial) <= (first if not lines else avail):
            line = trial
        else:
            lines.append(line)
            line = word
    lines.append(line)
    plt.close(probe)
    return head, lines, width(head + " ", bold=True) / width_in


def new_figure(width, height, caption, number, **subplot_kw):
    """Figure with a caption band at the bottom: rule line, then the text."""
    head, lines, head_frac = _wrap(caption, number, width)
    cap_in = 0.17 + 0.165 * len(lines)
    total = height + cap_in
    fig, ax = plt.subplots(figsize=(width, total), layout="constrained", **subplot_kw)
    frac = cap_in / total
    fig.get_layout_engine().set(rect=(0, frac, 1, 1 - frac))
    fig.add_artist(Line2D([0.012, 0.988], [frac - 0.012 / total] * 2,
                          color=RULE, lw=0.8, transform=fig.transFigure))
    y0 = frac - 0.075 / total
    kw = dict(va="top", fontsize=CAPTION_SIZE, color=INK)
    fig.text(0.012, y0, head, fontweight="bold", **kw)
    fig.text(0.012 + head_frac, y0, lines[0], **kw)
    if len(lines) > 1:
        fig.text(0.012, y0 - 0.165 / total, "\n".join(lines[1:]), linespacing=1.45, **kw)
    return fig, ax


def panel_label(a, letter, text=""):
    a.text(0.0, 1.03, letter, transform=a.transAxes, fontweight="bold", fontsize=10,
           color=INK, va="bottom")
    if text:
        a.text(0.055, 1.03, text, transform=a.transAxes, fontsize=9, color=GRAPHITE,
               va="bottom")


def end_label(a, x, y, text, color, dx=4, dy=0, **kw):
    a.annotate(text, (x, y), xytext=(dx, dy), textcoords="offset points", color=color,
               va="center", fontsize=8.5, **kw)


def save(fig, filename):
    fig.savefig(OUT / filename, dpi=200)
    plt.close(fig)


def _span(t_start, t_end):
    pad = 0.03 * (t_end - t_start)
    return t_start - pad, t_end + pad


# --- Figures 4-6: entry profiles -------------------------------------------

def entry_figure(entry, filename, number, caption, t_end, late, g_flight, g_limit=None,
                 q_limit=None, q_max=None, g_max=None, h_range=(0, 140), t_start=0.0,
                 rad_upper=None, design=None, legend_at=(1.0, 1.0), late_dy=0):
    """(a) altitude vs airspeed, (b) heating vs time, (c) deceleration vs time.

    late: (time s, flight altitude m) marked on (a).
    g_flight: flight peak deceleration, Earth g.
    """
    tr = simulate(**entry, t_end=t_end, dt_out=0.1)
    g = tr.decel / G0
    fig, ax = new_figure(10.5, 3.3, caption, number, ncols=3)
    xl = _span(t_start, t_end)

    a = ax[0]
    a.plot(tr.v / 1e3, tr.h / 1e3, color=INK)
    t_late, h_late = late
    k = np.argmin(abs(tr.t - t_late))
    a.plot(tr.v[k] / 1e3, h_late / 1e3, "o", ms=7, mfc=BG, mec=INK, mew=1.5)
    end_label(a, tr.v[k] / 1e3, h_late / 1e3, f"flight, {t_late:.0f} s", GRAPHITE, dx=7,
              dy=late_dy)
    a.set(xlabel="Airspeed (km/s)", ylabel="Altitude (km)", ylim=h_range,
          xlim=(0, 1.04 * tr.v.max() / 1e3))
    panel_label(a, "a")

    a = ax[1]
    total = np.where(np.isnan(tr.q_rad), np.nan, tr.q_total)
    a.plot(tr.t, total / 1e4, color=INK, lw=1.4, label=f"total {np.nanmax(total) / 1e4:.0f}")
    a.plot(tr.t, tr.q_conv / 1e4, color=SLATE, ls=DASH,
           label=f"convective {tr.peak_q_conv / 1e4:.0f}")
    a.plot(tr.t, tr.q_rad / 1e4, color=RED, label=f"radiative {tr.peak_q_rad / 1e4:.0f}")
    if rad_upper is not None:
        up = simulate(**{**entry, "body": replace(entry["body"], radiative=rad_upper)},
                      t_end=t_end, dt_out=0.1)
        a.plot(up.t, up.q_rad / 1e4, color=RED, ls=DASH, lw=1.4,
               label=f"older CN {up.peak_q_rad / 1e4:.0f}")
    if design is not None:
        a.plot(design["t"], design["q_conv"], "o", ms=6, mfc=BG, mec=SLATE, mew=1.4)
        a.plot(design["t"], design["q_rad_nonboltzmann"], "o", ms=6, mfc=BG, mec=RED, mew=1.4)
        a.plot(design["t"], design["q_rad_boltzmann"], "s", ms=5.5, mfc=BG, mec=RED, mew=1.4)
    a.legend(loc="upper right", bbox_to_anchor=legend_at, fontsize=8, handlelength=1.8,
             labelspacing=0.3, borderaxespad=0.2, handletextpad=0.5,
             title="peak, W/cm²", title_fontsize=8, alignment="left")
    if q_limit is not None:
        a.axhline(q_limit, color=GRAPHITE, lw=1, ls=DASH)
        end_label(a, xl[0], q_limit, f"limit {q_limit:.0f}", GRAPHITE, dx=2, dy=6)
    top = q_max or 1.15 * np.nanmax(total) / 1e4
    a.set(xlabel="Time from entry (s)", ylabel="Heating at the nose (W/cm²)", xlim=xl,
          ylim=(0, top))
    panel_label(a, "b")

    a = ax[2]
    a.plot(tr.t, g, color=INK)
    j = g.argmax()
    a.plot(tr.t[j], g_flight, "o", ms=7, mfc=BG, mec=INK, mew=1.5)
    g_text = f"{g_flight:.1f}".rstrip("0").rstrip(".")
    end_label(a, tr.t[j], g_flight, f"flight {g_text}", GRAPHITE, dx=8, dy=3)
    if g_limit is not None:
        a.axhline(g_limit, color=GRAPHITE, lw=1, ls=DASH)
        end_label(a, xl[1], g_limit, f"limit {g_limit:.0f}", GRAPHITE, dx=-2, dy=6, ha="right")
    a.set(xlabel="Time from entry (s)", ylabel="Deceleration (Earth g)", xlim=xl,
          ylim=(0, g_max or 1.25 * max(g.max(), g_limit or 0)))
    panel_label(a, "c")

    save(fig, filename)
    return tr


# --- Figures 7-8: corridor maps --------------------------------------------

def corridor_figure(entry, inertial, convert, filename, number, caption, speeds, angles_deg,
                    levels, point_label, g_limit=None, q_limit=None, g_label_at=None,
                    q_label_at=None):
    rot = dict(latitude=entry["latitude"], azimuth=entry["azimuth"])
    m = entry_corridor(entry["body"], entry["vehicle"], entry["h0"], speeds,
                       np.radians(angles_deg), convert=convert, **rot)
    edge = np.degrees(skip_out_boundary(
        entry["body"], entry["vehicle"], entry["h0"], speeds,
        np.radians(angles_deg[0]), np.radians(angles_deg[-1]), convert=convert, **rot))
    q, gload, V = m.peak_q_total / 1e4, m.peak_g, speeds / 1e3
    top = angles_deg[0]

    fig, a = new_figure(7.4, 4.6, caption, number)
    cmap = ListedColormap(HEAT_RAMP)
    cs = a.contourf(V, angles_deg, q, levels=levels, cmap=cmap,
                    norm=BoundaryNorm(levels, cmap.N), extend="max")
    a.contour(V, angles_deg, q, levels=levels[1:], colors=BG, linewidths=0.7)
    cb = fig.colorbar(cs, ax=a, pad=0.02, ticks=levels, aspect=30)
    cb.set_label("Peak heating at the nose (W/cm²)", color=GRAPHITE)
    cb.outline.set_visible(False)
    cb.ax.tick_params(length=0)

    if g_limit is not None:
        cg = a.contour(V, angles_deg, gload, levels=[g_limit], colors=[INK], linewidths=1.6,
                       linestyles=[DASH])
        a.clabel(cg, fmt=f"{g_limit:.0f} g", inline=True, fontsize=8.5, manual=[g_label_at])
    if q_limit is not None:
        cq = a.contour(V, angles_deg, q, levels=[q_limit], colors=[INK], linewidths=1.6,
                       linestyles=[DASH])
        a.clabel(cq, fmt=f"{q_limit:.0f} W/cm²", inline=True, fontsize=8.5,
                 manual=[q_label_at])

    if not np.all(np.isnan(edge)):
        a.fill_between(V, np.where(np.isnan(edge), top, edge), top, color=BG, zorder=3,
                       lw=0)
        a.plot(V, edge, color=INK, lw=1.2, zorder=4)
        a.text(V[-1] - 0.02 * (V[-1] - V[0]), top - 0.04 * (top - angles_deg[-1]),
               "skips out", ha="right", va="center", color=GRAPHITE, zorder=5)

    v_pt, g_pt = inertial["v"] / 1e3, np.degrees(inertial["gamma"])
    a.plot(v_pt, g_pt, "o", ms=8, mfc=BG, mec=INK, mew=1.6, zorder=6)
    a.annotate(point_label, (v_pt, g_pt), xytext=(-9, 0), textcoords="offset points",
               ha="right", va="center", color=INK, fontsize=9, zorder=6)
    a.set(xlabel="Entry speed, space-fixed (km/s)",
          ylabel="Entry flight-path angle, space-fixed (deg)",
          xlim=(V[0], V[-1]), ylim=(angles_deg[-1], top))
    a.grid(False)
    save(fig, filename)
    return m


# --- Figure 1: three worlds ------------------------------------------------

def three_worlds(runs, filename, number, caption):
    """runs: (world, mission, speed km/s, trajectory, optional upper-radiative trajectory)."""
    fig, ax = new_figure(8.6, 2.9, caption, number, ncols=len(runs))
    for a, letter, (world, mission, speed, tr, upper) in zip(ax, "abc", runs):
        conv, rad = tr.peak_q_conv / 1e4, tr.peak_q_rad / 1e4
        a.set_axisbelow(True)
        a.bar(0, conv, width=0.6, color=BG, edgecolor=SLATE, lw=1.6, hatch="////")
        a.bar(1, rad, width=0.6, color=RED, edgecolor=RED, lw=1.6)
        top = max(conv, rad)
        if upper is not None:
            hi = upper.peak_q_rad / 1e4
            a.plot([1, 1], [rad, hi], color=INK, lw=1)
            a.plot([0.88, 1.12], [hi, hi], color=INK, lw=1)
            top = max(top, hi)
            a.text(1.0, hi + 0.03 * top, f"{hi:.0f}", ha="center", va="bottom", color=INK,
                   fontsize=8.5)
        for x, val in ((0, conv), (1, rad)):
            y = val + 0.03 * top
            if upper is not None and x == 1:
                a.text(1.36, val, f"{val:.0f}", ha="left", va="center", color=INK, fontsize=8.5)
            else:
                a.text(x, y, f"{val:.0f}", ha="center", va="bottom", color=INK, fontsize=8.5)
        a.set_xticks([0, 1], ["convective", "radiative"])
        a.set_xlim(-0.6, 1.75)
        a.set_ylim(0, 1.18 * top)
        a.grid(axis="x", visible=False)
        a.tick_params(axis="x", length=0, labelcolor=INK)
        panel_label(a, letter, f"{world} · {mission} · {speed} km/s")
    ax[0].set_ylabel("Peak heating at the nose (W/cm²)")
    save(fig, filename)


# --- Figure 2: size effect -------------------------------------------------

def radiative_share(entry, nose_radius, t_end):
    vehicle = replace(entry["vehicle"], nose_radius=nose_radius)
    tr = simulate(**{**entry, "vehicle": vehicle}, t_end=t_end, dt_out=0.1)
    return tr.peak_q_rad / (tr.peak_q_rad + tr.peak_q_conv)


def size_figure(cases, titan, filename, number, caption):
    """cases: (world, color, entry, t_end). titan: (radius, low share, high share)."""
    fig, a = new_figure(7.4, 3.6, caption, number)
    a.axhline(50, color=GRAPHITE, lw=0.9, ls=DASH)
    crossings = {}
    for world, color, entry, t_end in cases:
        flown = entry["vehicle"].nose_radius
        lo, hi = entry["body"].radiative.valid["nose_radius"]
        radii = np.geomspace(lo, hi, 40)
        share = np.array([100 * radiative_share(entry, r, t_end) for r in radii])
        a.plot(radii, share, color=color)
        f = 100 * radiative_share(entry, flown, t_end)
        if flown < lo:
            a.plot([flown, lo], [f, share[0]], color=color, ls=(0, (1.5, 2)))
        a.plot(flown, f, "o", color=color, ms=6.5, mec=BG, mew=1.2)
        cross = np.interp(50, share, radii)
        crossings[world] = cross
        a.plot(cross, 50, "o", ms=6.5, mfc=BG, mec=color, mew=1.6)
        a.annotate(f"{cross:.1f} m", (cross, 50), xytext=(7, -6), textcoords="offset points",
                   ha="left", va="top", color=color, fontsize=8.5)
        end_label(a, radii[-1], share[-1], world, color, dx=6)
    radius, low, high = titan
    a.errorbar(radius, (low + high) / 2, yerr=(high - low) / 2, fmt="none", ecolor=IRIS,
               elinewidth=2, capsize=5, capthick=2)
    end_label(a, radius, high, "Titan, Huygens\n" + f"{low:.0f}–{high:.0f}%", IRIS,
              dx=-9, dy=-8, ha="right")
    a.set_xscale("log")
    a.set_xticks([0.1, 0.3, 1, 3, 10, 30], ["0.1", "0.3", "1", "3", "10", "30"])
    a.minorticks_off()
    a.set(xlim=(0.1, 60), ylim=(0, 75), xlabel="Nose radius (m)",
          ylabel="Radiative share of peak heating (%)")
    save(fig, filename)
    return crossings


# --- Figure 3: validation --------------------------------------------------

def validation_figure(cases, filename, number, caption, n=400):
    """cases: (label, entry, check time s, flight peak decel m/s^2, flight altitude m)."""
    fig, ax = new_figure(8.6, 2.5, caption, number, ncols=2, sharey=True)
    rows = range(len(cases))[::-1]
    for y, (label, entry, t_check, decel, altitude) in zip(rows, cases):
        peak, alt = disperse(entry, t_check, n=n)
        for a, err, unit in ((ax[0], 100 * (peak / decel - 1), "%"),
                             (ax[1], (alt - altitude) / 1e3, " km")):
            lo, hi = np.percentile(err, [2.5, 97.5])
            a.plot([lo, hi], [y, y], color=RULE, lw=6, solid_capstyle="butt")
            a.plot(err[0], y, "o", color=IRIS, ms=7, mec=BG, mew=1.2)
            a.annotate(f"{err[0]:+.1f}{unit}".replace("-", "−"), (err[0], y), xytext=(0, 9),
                       textcoords="offset points", ha="center", color=INK, fontsize=8.5,
                       bbox=dict(boxstyle="square,pad=0.1", fc=BG, ec="none"), zorder=5)
    for a, letter, xlabel in ((ax[0], "a", "Peak deceleration, model − flight (%)"),
                              (ax[1], "b", "Late-entry altitude, model − flight (km)")):
        a.axvline(0, color=INK, lw=1)
        a.set(xlim=(-8, 8), ylim=(-0.6, len(cases) - 0.35), xlabel=xlabel)
        a.grid(axis="y", visible=False)
        a.tick_params(axis="y", length=0, labelcolor=INK)
        panel_label(a, letter)
    ax[0].set_yticks(list(rows), [c[0] for c in cases])
    save(fig, filename)


def report(name, tr):
    print(f"{name}: peak {tr.peak_q_conv / 1e4:.0f} convective + {tr.peak_q_rad / 1e4:.0f} "
          f"radiative, {tr.peak_q_total / 1e4:.0f} total W/cm2; peak {tr.peak_decel_g:.1f} g")


def share_pct(conv_tr, rad_tr):
    return 100 * rad_tr.peak_q_rad / (rad_tr.peak_q_rad + conv_tr.peak_q_conv)


def main(corridors=True):
    OUT.mkdir(parents=True, exist_ok=True)
    SF, PF, HF = STARDUST_FLIGHT, PATHFINDER_FLIGHT, HUYGENS_FLIGHT

    sd = entry_figure(
        STARDUST_ENTRY, "stardust_entry.png", 4,
        "Stardust at Earth, entering at 12.9 km/s and −8.2°. (a) Altitude against airspeed; "
        "open circle: flight altitude at drogue deploy. (b) Heating at the nose; dashed: "
        "1200 W/cm² design limit. (c) Deceleration; open circle: flight peak; dashed: "
        "40 g design limit.",
        t_end=140, late=(SF["drogue_time"], SF["drogue_altitude"]), g_flight=SF["peak_g"],
        g_limit=40.0, q_limit=SF["peak_q_design"], q_max=1300, g_max=44, legend_at=(1.0, 0.88))
    report("Stardust", sd)

    pf = entry_figure(
        PATHFINDER_ENTRY, "pathfinder_entry.png", 5,
        "Mars Pathfinder, entering at 7.26 km/s and −14.1°. (a) Altitude against airspeed; "
        "open circle: flight altitude at parachute deploy (approximate). (b) Heating at the "
        "nose. (c) Deceleration; open circle: flight peak, reported as 16 g.",
        t_end=PF["chute_time"], late=(PF["chute_time"], PF["chute_altitude"]),
        g_flight=PF["peak_g"], q_max=135, g_max=20)
    report("Pathfinder", pf)

    hy = entry_figure(
        HUYGENS_ENTRY, "huygens_entry.png", 6,
        "Huygens at Titan, entering at 6.03 km/s and −65.5°. (a) Altitude against airspeed; "
        "open circle: last record of the measured entry profile. (b) Heating at the nose, "
        "with NASA design-study points (blue: convective; red circles and squares: "
        "radiative, newer and older CN models). Radiative lines cover only the speed range of "
        "those points. (c) Deceleration; open circle: flight peak, read from a published plot.",
        t_end=HF["late_time"], late=(HF["late_time"], HF["late_altitude"]),
        g_flight=HF["peak_decel"] / G0, q_max=95, g_max=14, h_range=(0, 1300),
        t_start=140.0, rad_upper=TITAN_HUYGENS_BOLTZMANN, design=HUYGENS_DESIGN_POINTS,
        legend_at=(1.0, 0.78), late_dy=12)
    report("Huygens", hy)
    hy_up = simulate(**{**HUYGENS_ENTRY, "body": replace(TITAN, radiative=TITAN_HUYGENS_BOLTZMANN)},
                     t_end=HF["late_time"], dt_out=0.1)

    three_worlds(
        [("Earth", "Stardust", 12.9, sd, None), ("Mars", "Pathfinder", 7.3, pf, None),
         ("Titan", "Huygens", 6.0, hy, hy_up)],
        "three_worlds.png", 1,
        "Peak heating at the nose as flown: convective (hatched) and radiative (filled), "
        f"W/cm². Radiative share: Earth {share_pct(sd, sd):.0f}%, Mars {share_pct(pf, pf):.0f}%, "
        f"Titan {share_pct(hy, hy):.0f}–{share_pct(hy, hy_up):.0f}%. Titan's whisker runs to "
        "the older CN radiation model.")

    size_figure(
        [("Earth", SLATE, STARDUST_ENTRY, 140), ("Mars", BRASS, PATHFINDER_ENTRY, PF["chute_time"])],
        (HUYGENS_ENTRY["vehicle"].nose_radius, share_pct(hy, hy), share_pct(hy, hy_up)),
        "size_effect.png", 2,
        "Radiative share of peak heating against nose radius, with each trajectory fixed: "
        "Stardust's at Earth, Pathfinder's at Mars. Solid: radii the radiative fit covers; "
        "dotted: extension to the flown capsule (filled dot). Above the dashed line radiation "
        "exceeds convection. Titan: one point at Huygens' 1.25 m nose, a range because two CN "
        "radiation models disagree; no curve, since Titan heating data exists for that size only.")

    s = SIGMA
    validation_figure(
        [("Stardust", STARDUST_ENTRY, SF["drogue_time"], SF["peak_g"] * G0, SF["drogue_altitude"]),
         ("Pathfinder", PATHFINDER_ENTRY, PF["chute_time"], PF["peak_g"] * G0,
          PF["chute_altitude"]),
         ("Huygens", HUYGENS_ENTRY, HF["late_time"], HF["peak_decel"], HF["late_altitude"])],
        "validation.png", 3,
        "Model minus flight. Dot: nominal run. Bar: central 95% of 400 runs with entry angle, "
        f"entry speed, drag coefficient and density perturbed (σ = {np.degrees(s['gamma']):.2f}°, "
        f"{100 * s['speed']:.1f}%, {100 * s['drag']:.0f}%, {100 * s['density']:.0f}%). "
        "(a) Peak deceleration. (b) Altitude late in entry: drogue deploy (Stardust), "
        "parachute deploy (Pathfinder), last measured entry record (Huygens).")

    if corridors:
        corridor_figure(
            STARDUST_ENTRY, STARDUST_INERTIAL, stardust_convert, "stardust_corridor.png", 7,
            "Peak heating for the Stardust capsule over entry speed and angle. Above the top "
            "curve it skips back out. Dashed: limits; below the 40 g curve it exceeds the "
            "deceleration limit, right of the 1200 W/cm² curve the heating limit. Open circle: "
            "Stardust's entry.",
            speeds=np.linspace(8e3, 16e3, 41), angles_deg=np.linspace(-4.0, -16.0, 49),
            levels=[0, 200, 400, 600, 800, 1200, 1600, 2400], point_label="Stardust",
            g_limit=40.0, q_limit=1200.0, g_label_at=(9.6, -12.2), q_label_at=(12.4, -13.0))
        corridor_figure(
            PATHFINDER_ENTRY, PATHFINDER_INERTIAL, pathfinder_convert, "pathfinder_corridor.png",
            8,
            "Peak heating for the Pathfinder capsule at Mars over entry speed and angle. Above "
            "the top curve it skips back out. Open circle: Pathfinder's entry.",
            speeds=np.linspace(4e3, 9e3, 41), angles_deg=np.linspace(-8.0, -24.0, 49),
            levels=[0, 25, 50, 75, 100, 150, 200, 300], point_label="Pathfinder")
    print("Figures in", OUT)


if __name__ == "__main__":
    main()
