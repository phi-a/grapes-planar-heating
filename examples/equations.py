"""Figure 1: the two heating models, typeset by LaTeX with color-coded inputs.

Needs LaTeX with xcolor, plus dvipng. Run: python examples/equations.py
"""
import shutil
import subprocess
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg

from showcase import BRASS, GRAPHITE, INK, IRIS, SLATE, new_figure, save

INPUTS = [  # (symbol, role, color)
    (r"$\rho$", "atmosphere", BRASS),
    (r"$V$", "trajectory", IRIS),
    (r"$R_n$", "vehicle", SLATE),
    (r"$k, C, a, b, f$", "gas chemistry", GRAPHITE),
]

PREAMBLE = r"""\documentclass[varwidth]{standalone}
\usepackage{amsmath}
\usepackage{xcolor}
\definecolor{ink}{HTML}{%s}
\definecolor{atm}{HTML}{%s}
\definecolor{trj}{HTML}{%s}
\definecolor{veh}{HTML}{%s}
\definecolor{gas}{HTML}{%s}
\newcommand{\rh}{{\color{atm}\rho}}
\newcommand{\V}{{\color{trj}V}}
\newcommand{\Rn}[1][]{{\color{veh}R}_{\color{veh}n}^{#1}}
\newcommand{\g}[1]{{\color{gas}#1}}
\begin{document}\color{ink}\Large
""" % tuple(c.lstrip("#") for c in (INK, BRASS, IRIS, SLATE, GRAPHITE))

EQUATIONS = [
    ("Convective\nSutton–Graves",
     r"$\displaystyle q_c = \g{k}\,\sqrt{\frac{\rh}{\Rn}}\;\V^{3}$"),
    ("Radiative\nTauber–Sutton",
     r"$\displaystyle q_r = \g{C}\;\Rn[\g{a}]\;\rh^{\g{b}}\;\g{f}(\V)$"),
    ("Their ratio",
     r"$\displaystyle \frac{q_r}{q_c} \;\propto\; \Rn[\g{a}+\frac{1}{2}]\;"
     r"\rh^{\g{b}-\frac{1}{2}}\;\frac{\g{f}(\V)}{\V^{3}}$"),
]

DPI = 1000  # typesetting resolution; images are placed at their natural size


def render(tex, workdir, name):
    """Typeset one formula with latex + dvipng on a white background."""
    src = Path(workdir) / f"{name}.tex"
    src.write_text(PREAMBLE + tex + "\n\\end{document}\n", encoding="utf-8")
    run = dict(cwd=workdir, check=True, capture_output=True)
    subprocess.run(["latex", "-interaction=nonstopmode", src.name], **run)
    subprocess.run(["dvipng", "-D", str(DPI), "-T", "tight", "-bg", "White", "--truecolor",
                    "-o", f"{name}.png", f"{name}.dvi"], **run)
    return mpimg.imread(Path(workdir) / f"{name}.png")


def main():
    if not (shutil.which("latex") and shutil.which("dvipng")):
        raise SystemExit("needs latex and dvipng on PATH")
    caption = (
        "Stagnation-point heating models used here, SI units; colors mark where each input "
        "comes from. Convection falls with nose radius; radiation rises with it, with density "
        "(b ≈ 1.2) and steeply with speed: at Earth f(V) grows about as V⁷ between 11 and "
        "16 km/s. The ratio shows why radiation takes over for larger, denser and faster "
        "entries. Titan uses a separate fit to published points.")
    fig, ax = new_figure(7.4, 3.0, caption, 1)
    ax.set_axis_off()
    with tempfile.TemporaryDirectory() as tmp:
        images = [render(tex, tmp, f"eq{i}") for i, (_, tex) in enumerate(EQUATIONS)]

    fig.canvas.draw()
    box = ax.get_position()
    W, H = fig.get_figwidth(), fig.get_figheight()
    for (label, _), img, y in zip(EQUATIONS, images, (0.88, 0.6, 0.3)):
        h, w = img.shape[0] / DPI / H, img.shape[1] / DPI / W
        yc = box.y0 + y * box.height
        a = fig.add_axes([box.x0 + 0.25 * box.width, yc - h / 2, w, h])
        a.imshow(img, interpolation="lanczos")
        a.set_axis_off()
        fig.text(box.x0 + 0.01, yc, label, ha="left", va="center", fontsize=9.5, color=INK,
                 linespacing=1.3)

    x, y = box.x0 + 0.01, box.y0 + 0.035
    for symbol, role, color in INPUTS:
        t = fig.text(x, y, f"{symbol}  {role}", va="center", fontsize=9, color=color)
        fig.canvas.draw()
        x = t.get_window_extent().transformed(fig.transFigure.inverted()).x1 + 0.04
    save(fig, "equations.png")


if __name__ == "__main__":
    main()
