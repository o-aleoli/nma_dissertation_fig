#!/usr/bin/env python3
import numpy as np
from matplotlib.backends.backend_pgf import PdfPages
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
import os

input_files_Nb = [
    'nn_total',
    'Nb_d_separated_internal'
]
input_files_Ta = [
    'nt_total',
    'Ta_d_separated_internal'
]
output_file = 'nt_nn_internal_pdos'
xlim = (1, 7)
ylim = (-10, 10)

data_Nb = []
for i in input_files_Nb:
    data_Nb.append(np.load(i+'.npz'))

data_Ta = []
for i in input_files_Ta:
    data_Ta.append(np.load(i+'.npz'))

width = 210 - 30 - 20
height = 85

"""
INPUT_DATA is a pickled, compressed file with the following numpy arrays:
    * pdos, an np.array, with dimensions (ispin, numpoints), where ispin=2 and numpoints=301 (typically)
    * eigen_energy, an np.array, with dimensions (numpoints)
    * eigen_adjust, a float
    * e_fermi, a float
"""
with PdfPages(output_file+'.pdf') as output:
    plt.rcParams.update({
        'font.family': 'serif',
        'font.size': 11.0,
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
        nrows=2,
        ncols=1,
        constrained_layout=True
    )

    d_colors = [
        "#0d8200",#dxy
        "#13bb00",#dyz
        "#ff00d9",#dz2
        "#18f500",#dxz
        "#860069"# dx2
    ]

    legend = [
        Line2D(
            [0], [0],
            color=d_colors[4],
            linewidth=1.0,
            label=r'M $d_{x^2 - y^2}$'
        ),
        Line2D(
            [0], [0],
            color=d_colors[2],
            linewidth=1.0,
            label=r'M $d_{3z^2 - r^2}$'
        ),
        Line2D(
            [0], [0],
            color=d_colors[0],
            linewidth=1.0,
            label=r'M $d_{xy}$'
        ),
        Line2D(
            [0], [0],
            color=d_colors[1],
            linewidth=1.0,
            label=r'M $d_{yz}$'
        ),
        Line2D(
            [0], [0],
            color=d_colors[3],
            linewidth=1.0,
            label=r'M $d_{xz}$'
        ),
        Line2D(
            [0], [0],
            color='black',
            linewidth=1.0,
            label='Total'
        ),
        Line2D(
            [0], [0],
            color='black',
            linewidth=2.0,
            linestyle='dashed',
            label=r'$\varepsilon_F$'
        )
    ]

###############################################################################
    alignment_Ta = 1.936701

    x_axis = data_Ta[0]['eigen_energy'] + alignment_Ta

    axs[0].plot(
        x_axis,  data_Ta[0]['pdos'][0, ...],
        x_axis, -data_Ta[0]['pdos'][0, ...],
        c='black',
        linewidth=1.0
    )

    for orb in range(len(d_colors)):
        axs[0].plot(
            x_axis,  data_Ta[1]['pdos'][0, ...,4 + orb],
            x_axis, -data_Ta[1]['pdos'][1, ...,4 + orb],
            linewidth=1.0,
            c=d_colors[orb]
        )

    axs[0].plot(
        [data_Ta[0]['e_fermi'] + alignment_Ta, data_Ta[0]['e_fermi'] + alignment_Ta],
        [np.amax(data_Ta[0]['pdos'][0, ...]), -np.amax(data_Ta[0]['pdos'][0, ...])],
        linestyle='dashed',
        linewidth=1.0,
        c='black'
    )

###############################################################################
    alignment_Nb = 2.231016

    x_axis = data_Nb[0]['eigen_energy'] + alignment_Nb

    axs[1].plot(
        x_axis,  data_Nb[0]['pdos'][0, ...],
        x_axis, -data_Nb[0]['pdos'][1, ...],
        c='black',
        linewidth=1.0
    )

    for orb in range(len(d_colors)):
        axs[1].plot(
            x_axis,  data_Nb[1]['pdos'][0, ...,4 + orb],
            x_axis, -data_Nb[1]['pdos'][1, ...,4 + orb],
            linewidth=1.0,
            c=d_colors[orb]
        )

    axs[1].plot(
        [data_Nb[0]['e_fermi'] + alignment_Nb, data_Nb[0]['e_fermi'] + alignment_Nb],
        [np.amax(data_Nb[0]['pdos'][0, ...]), -np.amax(data_Nb[0]['pdos'][0, ...])],
        linestyle='dashed',
        linewidth=1.0,
        c='black'
    )

###############################################################################

    axs[0].set(
        xlim=xlim,
        ylim=(-0.5, 0.5),
        # yticks=list(ylim),
        yticks=[],
        xticks=[],
        ylabel=r"pDOS (eV${}^{-1}$)"
    )

    axs[1].set(
        xlim=xlim,
        ylim=ylim,
        yticks=[],
        xlabel=r"Energy (eV)",
        ylabel=r"pDOS (eV${}^{-1}$)"
    )

    axs[0].text(0.05, 0.8, 'a', transform=axs[0].transAxes)
    axs[1].text(0.05, 0.8, 'b', transform=axs[1].transAxes)

    axs[1].legend(
        handles=legend,
        loc='lower center',
        bbox_to_anchor=(0.485, -0.9),
        ncol=6,
        fontsize='small'
    )
    
    output.savefig()
