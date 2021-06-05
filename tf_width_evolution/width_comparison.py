#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pgf import PdfPages
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

data = [
   pd.read_fwf("1cell.dat"),
   pd.read_fwf("3cell.dat"),
   pd.read_fwf("5cell.dat")
]
color = [
    "black",
    "blue",
    "green"
]
label = [
    "One-cell",
    "Three-cell",
    "Five-cell"
]
surface = pd.read_fwf("surface.dat")

width = 210 - 30 - 20 - 25
height = 80
output_file = "width_comparison"

with PdfPages(output_file + '.pdf') as output:
    plt.rcParams.update({
    "font.family": "serif",
    "font.size": "11",
    "figure.figsize": [width/25.4, height/25.4],
    "savefig.dpi": 1200,
    "text.usetex": True,
    "pgf.rcfonts": False,
    "pgf.texsystem": "pdflatex",
    "pgf.preamble": r"\usepackage{mathpazo,eulervm}\usepackage[utf8x]{inputenc}\usepackage[separate-uncertainty = true]{siunitx}",
    "legend.frameon": False,
    "legend.fancybox": False,
    "legend.fontsize": "small"
    })

    fig, axs = plt.subplots(
        ncols=2,
        sharey=True,
        constrained_layout=True
    )

    for slab in range(len(data)):
        axs[0].scatter(
            data[slab]["bond"],
            data[slab]["layer"],
            marker = 'o',
            color = color[slab]
        )
        axs[1].scatter(
            data[slab]["baur_atom"],
            data[slab]["layer"],
            marker = "s",
            color = color[slab],
            label = label[slab]
        )

    axs[0].scatter(
            surface["bond"],
            surface["layer"],
            marker = "x",
            color = "red",
            label = "Inf. surface"
    )

    axs[1].scatter(
            surface["baur_atom"],
            surface["layer"],
            marker = "+",
            color = "red",
            label = "Inf. surface"
    )

    axs[0].set(
        xlim = [1.90, 2.10],
        ylabel = "Layers",
        yticks = [],
        xlabel = r"$d_{\text{Nb--O}}$ (\si{\angstrom})",
    )

    axs[0].xaxis.set_minor_locator(AutoMinorLocator())

    axs[1].set(
        xlim = [-0.5, 7.5],
        xlabel = r"$\Delta_d$ (\si{\percent})",
    )

    axs[1].xaxis.set_minor_locator(AutoMinorLocator())

    axs[1].legend()

    output.savefig()
