#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import scipy.optimize as spop
from matplotlib.backends.backend_pgf import PdfPages


def scaling_law(x, a0, a1, a3):
    return a0 + a1*x + a3*x**3


# Graph inputs
width = 210.0 - 20.0 - 30.0
height = 120.0
output_file = "fit_formation_energy_E"

data_Na0 = np.array([
    [0.15956992, 4.05355742],
    [0.0797855, 4.27098162],
    [0.05319029, 4.66385604]
])

data_O0 = np.array([
    [0.15956992,  1.14283325],
    [0.0797855, -0.61898431],
    [0.05319029, -0.21404431]
])

data_Na1 = np.array([
    [0.15956992, 5.03448],
    [0.0797855, 5.54564],
    [0.05319029, 5.96716]
])

data_O2 = np.array([
    [0.15956992, -2.40674],
    [0.0797855, -3.05966],
    [0.05319029, -2.79944]
])

with PdfPages(output_file + ".pdf") as output:
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 11,
        "figure.figsize": [width/25.4, height/25.4],
        "savefig.dpi": 1200,
        "text.usetex": True,
        "pgf.rcfonts": False,
        "pgf.texsystem": "pdflatex",
        "pgf.preamble": r"\usepackage{mathpazo,eulervm,xcolor}\usepackage[utf8x]{inputenc}\usepackage[separate-uncertainty = true]{siunitx}",
        "legend.frameon": False,
        "legend.fancybox": False,
        "legend.fontsize": "small"
    })
    fig, axs = plt.subplots(
        nrows = 2,
        ncols = 2,
        constrained_layout = True
    )

    axs[0, 0].plot(
        np.linspace(0, data_Na0[0, 0], num = 100),
        scaling_law(
            np.linspace(0, data_Na0[0, 0], num = 100),
            5.668059108160884,
            -19.974657964061887,
            387.11056012780386
        ),
        color = 'black',
        linestyle = 'dashed',
        linewidth = 1.0
    )
    axs[0, 0].plot(
        data_Na0[..., 0],
        data_Na0[..., 1],
        "ko",
    )
    axs[0, 1].plot(
        np.linspace(0, data_O0[0, 0], num = 100),
        scaling_law(
            np.linspace(0, data_O0[0, 0], num = 100),
            1.2723558994657198,
            -31.33665452096673,
            1198.8172575121782
        ),
        color = 'black',
        linestyle = 'dashed',
        linewidth = 1.0
    )
    axs[0, 1].plot(
        data_O0[..., 0],
        data_O0[..., 1],
        "ko",
    )
    axs[1, 0].plot(
        np.linspace(0, data_Na1[0, 0], num = 100),
        scaling_law(
            np.linspace(0, data_Na1[0, 0], num = 100),
            6.981425514192862,
            -19.92705857505196,
            303.42009724162085
        ),
        color = 'black',
        linestyle = 'dashed',
        linewidth = 1.0
    )
    axs[1, 0].plot(
        data_Na1[..., 0],
        data_Na1[..., 1],
        "ko",
    )
    axs[1, 1].plot(
        np.linspace(0, data_O2[0, 0], num = 100),
        scaling_law(
            np.linspace(0, data_O2[0, 0], num=100),
            -1.9531817038097423,
            -17.54348987364265,
            577.3618337129745
        ),
        color = 'black',
        linestyle = 'dashed',
        linewidth = 1.0
    )
    axs[1, 1].plot(
        data_O2[..., 0],
        data_O2[..., 1],
        "ko",
    )
    labels = [
        Line2D(
            [0], [0],
            color = 'black',
            linewidth = 1.0,
            linestyle = 'dashed',
            label = r"Fitted Eq.~2.14"
        )
    ]
    axs[1, 1].legend(
        handles = labels,
        loc = 'lower center',
        bbox_to_anchor = (0.75, -0.5),
        ncol = 1
    )
    axs[0, 0].set(
        ylabel = r"$\Delta E_f$ @ E (\si{\electronvolt})",
        xlim = [0, 0.175],
        xticks = (0.00, 0.05, 0.08, 0.16),
        xticklabels = []
    )
    secondary_00 = axs[0, 0].secondary_xaxis("top")
    secondary_00.set_xlim([0, 0.175])
    secondary_00.set_xticks([0.00, 0.05, 0.08, 0.16])
    secondary_00.set_xticklabels([r"$\infty$", r"$3\!\times\!3\!\times\!3$", r"$2\!\times\!2\!\times\!2$", r"$1\!\times\!1\!\times\!1$"], fontsize = 8)
    axs[0, 1].set(
        xlim = [0, 0.175],
        xticks = (0.00, 0.05, 0.08, 0.16),
        xticklabels = []
    )
    secondary_01 = axs[0, 1].secondary_xaxis("top")
    secondary_01.set_xlim([0, 0.175])
    secondary_01.set_xticks([0.00, 0.05, 0.08, 0.16])
    secondary_01.set_xticklabels([r"$\infty$", r"$3\!\times\!3\!\times\!3$", r"$2\!\times\!2\!\times\!2$", r"$1\!\times\!1\!\times\!1$"], fontsize = 8)
    axs[1, 0].set(
        xlim = [0, 0.175],
        xticks = (0.00, 0.05, 0.08, 0.16),
        xlabel = r"$L^{-1}$ (\si{\per\angstrom})",
        ylabel = r"$\Delta E_f$ @ E (\si{\electronvolt})"
    )
    axs[1, 1].set(
        xlim = [0, 0.175],
        xticks = (0.00, 0.05, 0.08, 0.16),
        xlabel = r"$L^{-1}$ (\si{\per\angstrom})",
    )

    axs[0, 0].annotate(
        r"\color{red}{\SI{5.67}{\electronvolt}}",
        xy = (0.0000, 5.668059108160884),
        xycoords = 'data',
        xytext = (0.15, 0.95),
        arrowprops = {
            'color': 'red',
            'width': 1.0,
            'headwidth': 4,
            'shrink': 0.05
        },
        horizontalalignment = 'left',
        verticalalignment = 'top',
        textcoords = 'axes fraction'
    )
    axs[1, 0].annotate(
        r"\color{red}{\SI{6.98}{\electronvolt}}",
        xy = (0.0000, 6.981425514192862),
        xycoords = 'data',
        xytext = (0.15, 0.95),
        arrowprops = {
            'color': 'red',
            'width': 1.0,
            'headwidth': 4,
            'shrink': 0.05
        },
        horizontalalignment = 'left',
        verticalalignment = 'top',
        textcoords = 'axes fraction'
    )

    axs[0, 1].annotate(
        r"\color{red}{\SI{1.27}{\electronvolt}}",
        xy = (0.0000, 1.2723558994657198),
        xycoords = 'data',
        xytext = (0.15, 0.95),
        arrowprops = {
            'color': 'red',
            'width': 1.0,
            'headwidth': 4,
            'shrink': 0.05
        },
        horizontalalignment = 'left',
        verticalalignment = 'top',
        textcoords = 'axes fraction'
    )
    axs[1, 1].annotate(
        r"\color{red}{\SI{-1.95}{\electronvolt}}",
        xy = (0.0000, -1.9531817038097423),
        xycoords = 'data',
        xytext = (0.15, 0.95),
        arrowprops = {
            'color': 'red',
            'width': 1.0,
            'headwidth': 4,
            'shrink': 0.05
        },
        horizontalalignment = 'left',
        verticalalignment = 'top',
        textcoords = 'axes fraction'
    )

    axs[0, 0].text(
        0.95,
        0.85,
        "a",
        transform = axs[0, 0].transAxes
    )
    axs[0, 1].text(
        0.95,
        0.85,
        "b",
        transform = axs[0, 1].transAxes
    )
    axs[1, 0].text(
        0.95,
        0.85,
        "c",
        transform = axs[1, 0].transAxes
    )
    axs[1, 1].text(
        0.95,
        0.85,
        "d",
        transform = axs[1, 1].transAxes
    )

    output.savefig()
