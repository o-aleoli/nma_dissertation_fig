#!/usr/bin/env python3
"""
input_files is a pickled, compressed file with the following numpy arrays:
    * pdos, an np.array, with dimensions (ispin, numpoints), where ispin=2 and numpoints=301 (typically)
    * eigen_energy, an np.array, with dimensions (numpoints)
    * eigen_adjust, a float
    * e_fermi, a float
"""
import numpy as np
from matplotlib.backends.backend_pgf import PdfPages
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

input_files = [
    "Nb_surface_separated_d",
    "Nb_subsurface_separated_d",
    "O_total_p",
    "total"
]
alignment = 2.231016
output_file = 'pdos_nn_thinfilm_surface_subsurface'
xlim = (-4, 12)
ylim = (-2.5, 2.5)
width = 210 - 30 - 20
height = 80

data = []
for i in input_files:
    data.append(np.load(i+'.npz'))

with PdfPages(output_file+'.pdf') as output:
    plt.rcParams.update({
        'font.family': 'serif',
        'font.size': 11,
        'figure.figsize': [width/25.4, height/25.4],
        'savefig.dpi': 1200,
        'text.usetex': True,
        'pgf.rcfonts': False,
        'pgf.texsystem': 'pdflatex',
        'pgf.preamble': r"\usepackage{mathpazo,eulervm}\usepackage[utf8x]{inputenc}\usepackage[separate-uncertainty = true]{siunitx}",
        'legend.frameon': False,
        'legend.fancybox': False
    })
    fig, axs = plt.subplots(
        nrows = 2,
        ncols = 1,
        constrained_layout = True
    )

    d_colors = [
        "#0d8200",#dxy
        "#13bb00",#dyz
        "#ff00d9",#dz2
        "#18f500",#dxz
        "#860069"# dx2
    ]

    LINES = [
        Line2D(
            [0], [0],
            color = d_colors[4],
            linewidth = 1.0,
            label = r'Nb $4d_{x^2 - y^2}$'
        ),
        Line2D(
            [0], [0],
            color = d_colors[2],
            linewidth = 1.0,
            label = r'Nb $4d_{3z^2 - r^2}$'
        ),
        Line2D(
            [0], [0],
            color = d_colors[0],
            linewidth = 1.0,
            label = r'Nb $4d_{xy}$'
        ),
        Line2D(
            [0], [0],
            color = d_colors[3],
            linewidth = 1.0,
            label = r'Nb $4d_{xz}$'
        ),
        Line2D(
            [0], [0],
            color = d_colors[1],
            linewidth = 1.0,
            label = r'Nb $4d_{zy}$'
        ),
        # Line2D(
        #     [0], [0],
        #     color = 'red',
        #     linewidth = 1.0,
        #     label = r'O $2p$'
        # ),
        Line2D(
            [0], [0],
            color = 'black',
            linewidth = 1.0,
            label = 'Total'
        ),
        Line2D(
            [0], [0],
            color = 'black',
            linewidth = 1.0,
            linestyle = 'dashed',
            label = r'$\varepsilon_F$'
        ),
    ]

###############################################################################

    x_axis = data[0]['eigen_energy'] + alignment 

    for orb in range(len(d_colors)):
        axs[0].plot(
            x_axis,  data[0]['pdos'][0, ...,4 + orb],
            x_axis, -data[0]['pdos'][1, ...,4 + orb],
            linewidth = 1.0,
            color = d_colors[orb]
        )

    # axs[0].plot(
    #     x_axis,  data[2]['pdos'][0, ...],
    #     x_axis, -data[2]['pdos'][1, ...],
    #     linewidth = 1.0,
    #     color = 'red'
    # )

    axs[0].plot(
        x_axis,  data[3]['pdos'][0, ...],
        x_axis, -data[3]['pdos'][1, ...],
        linewidth = 1.0,
        color = 'black'
    )

    axs[0].vlines(
        data[0]['e_fermi'] + alignment,
        ylim[0],
        ylim[1],
        linestyle = 'dashed',
        linewidth = 2.0,
        color = 'black'
    )

###############################################################################

    for orb in range(len(d_colors)):
        axs[1].plot(
            x_axis,  data[1]['pdos'][0, ...,4 + orb],
            x_axis, -data[1]['pdos'][1, ...,4 + orb],
            linewidth = 1.0,
            color = d_colors[orb]
        )

    # axs[1].plot(
    #     x_axis,  data[2]['pdos'][0, ...],
    #     x_axis, -data[2]['pdos'][1, ...],
    #     linewidth = 1.0,
    #     color = 'red'
    # )

    axs[1].plot(
        x_axis,  data[3]['pdos'][0, ...],
        x_axis, -data[3]['pdos'][1, ...],
        linewidth = 1.0,
        color = 'black'
    )

    axs[1].vlines(
        data[0]['e_fermi'] + alignment,
        ylim[0],
        ylim[1],
        linestyle = 'dashed',
        linewidth = 2.0,
        color = 'black'
    )

###############################################################################

    axs[1].legend(
        handles = LINES,
        loc = 'lower center',
        bbox_to_anchor = (0.5, -1.0),
        ncol = 4
    )

    axs[0].text(0.05, 0.79, 'a', transform = axs[0].transAxes)
    axs[1].text(0.05, 0.79, 'b', transform = axs[1].transAxes)

    axs[0].set(
        xlim = xlim,
        ylim = ylim,
        ylabel = r"pDOS (eV${}^{-1}$)",
        xticks = [],
        yticks = []
    )

    axs[1].set(
        xlim = xlim,
        ylim = ylim,
        ylabel = r"pDOS (eV${}^{-1}$)",
        xlabel = r"Energy (eV)",
        yticks = [],
    )

    output.savefig()
