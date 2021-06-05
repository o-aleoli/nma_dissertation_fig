#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pgf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

width = 210 - 30 - 20
height = 85
output_file = "width_pdos"
labels = ["One-cell", "Three-cell", "Five-cell", "Inf. surface"]
input_files_1cell = [
    "Nb_surface_separated_d_1cell",
    "Nb_subsurface_separated_d_1cell",
    "total_1cell"
]
input_files_3cell = [
    "Nb_surface_separated_d_3cell",
    "Nb_subsurface_separated_d_3cell",
    "total_3cell"
]
input_files_5cell = [
    "Nb_surface_separated_d_5cell",
    "Nb_subsurface_separated_d_5cell",
    "total_5cell"
]
input_files_surface = [
    "Nb_surface_separated_d_surface",
    "Nb_subsurface_separated_d_surface",
    "total_surface"
]
d_colors = [
    "#0d8200",#dxy
    "#13bb00",#dyz
    "#ff00d9",#dz2
    "#18f500",#dxz
    "#860069"# dx2
]
legends1 = [Line2D([0], [0], color = d_colors[4], linewidth = 1.0, label = r'Nb $4d_{x^2 - y^2}$')]
legends2 = [
    Line2D([0], [0], color = d_colors[2], linewidth = 1.0, label = r'Nb $4d_{3z^2 - r^2}$'),
    Line2D([0], [0], color = d_colors[0], linewidth = 1.0, label = r'Nb $4d_{xy}$')
]
legends3 = [
    Line2D([0], [0], color = d_colors[3], linewidth = 1.0, label = r'Nb $4d_{xz}$'),
    Line2D([0], [0], color = d_colors[1], linewidth = 1.0, label = r'Nb $4d_{zy}$')
]
legends4 = [
    Line2D([0], [0], color = 'black', linewidth = 2.0, linestyle = 'dashed', label = r'$\varepsilon_F$'),
    Line2D([0], [0], color = 'black', linewidth = 1.0, label = "Total")
]
xlim = (-4.5, 4.5)
ylim = (-2.5, 2.5)

data_1cell = []
for i in input_files_1cell:
    data_1cell.append(np.load(i + '.npz'))
data_3cell = []
for i in input_files_3cell:
    data_3cell.append(np.load(i + '.npz'))
data_5cell = []
for i in input_files_5cell:
    data_5cell.append(np.load(i + '.npz'))
data_surface = []
for i in input_files_surface:
    data_surface.append(np.load(i + '.npz'))

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

    _, axs = plt.subplots(
        nrows = 2,
        ncols = 4,
        constrained_layout = True,
        sharex = True,
        sharey = True
    )

    axs[0, 0].set(title = "One-cell", ylabel = r"pDOS (\si{\per\electronvolt})")
    axs[0, 1].set(title = "Three-cell")
    axs[0, 2].set(title = "Five-cell")
    axs[0, 3].set(title = "Inf. surface")
    
    axs[1, 0].set(xlabel = r"Energy (\si{\electronvolt})", ylabel = r"pDOS (\si{\per\electronvolt})")
    axs[1, 1].set(xlabel = r"Energy (\si{\electronvolt})")
    axs[1, 2].set(xlabel = r"Energy (\si{\electronvolt})")
    axs[1, 3].set(xlabel = r"Energy (\si{\electronvolt})")

    axs[1, 0].legend(handles=legends1, loc = "lower center", bbox_to_anchor = (0.5, -0.750))
    axs[1, 1].legend(handles=legends2, loc = "lower center", bbox_to_anchor = (0.5, -0.925))
    axs[1, 2].legend(handles=legends3, loc = "lower center", bbox_to_anchor = (0.5, -0.925))
    axs[1, 3].legend(handles=legends4, loc = "lower center", bbox_to_anchor = (0.5, -0.925))

    for i in range(2):
        x_axis = data_1cell[i]['eigen_energy']

        for orb in range(len(d_colors)):
            axs[i, 0].plot(
                x_axis,  data_1cell[i]['pdos'][0, :,4 + orb],
                x_axis, -data_1cell[i]['pdos'][1, :,4 + orb],
                linewidth=1.0,
                c=d_colors[orb]
            )

        axs[i, 0].plot(
            [data_1cell[i]['e_fermi'], data_1cell[i]['e_fermi']],
            [np.amax(data_1cell[0]['pdos'][0, :]), -np.amax(data_1cell[0]['pdos'][0, :])],
            linestyle = 'dashed',
            linewidth = 1.0,
            c = 'black'
        )

        axs[i, 0].plot(
            x_axis,  data_1cell[2]['pdos'][0, :],
            x_axis, -data_1cell[2]['pdos'][1, :],
            linewidth = 1.0,
            color = "black"
        )

        x_axis = data_3cell[i]['eigen_energy']

        for orb in range(len(d_colors)):
            axs[i, 1].plot(
                x_axis,  data_3cell[i]['pdos'][0, :,4 + orb],
                x_axis, -data_3cell[i]['pdos'][1, :,4 + orb],
                linewidth = 1.0,
                color = d_colors[orb]
            )

        axs[i, 1].plot(
            [data_3cell[i]['e_fermi'], data_3cell[i]['e_fermi']],
            [np.amax(data_3cell[0]['pdos'][0, :]), -np.amax(data_3cell[0]['pdos'][0, :])],
            linestyle = 'dashed',
            linewidth = 1.0,
            c = 'black'
        )

        axs[i, 1].plot(
            x_axis,  data_3cell[2]['pdos'][0, :],
            x_axis, -data_3cell[2]['pdos'][1, :],
            linewidth = 1.0,
            color = "black"
        )

        x_axis = data_5cell[i]['eigen_energy']

        for orb in range(len(d_colors)):
            axs[i, 2].plot(
                x_axis,  data_5cell[i]['pdos'][0, :,4 + orb],
                x_axis, -data_5cell[i]['pdos'][1, :,4 + orb],
                linewidth = 1.0,
                color = d_colors[orb]
            )

        axs[i, 2].plot(
            [data_5cell[i]['e_fermi'], data_5cell[i]['e_fermi']],
            [np.amax(data_5cell[i]['pdos'][0, :]), -np.amax(data_5cell[i]['pdos'][0, :])],
            linestyle = 'dashed',
            linewidth = 1.0,
            color = 'black'
        )

        axs[i, 2].plot(
            x_axis,  data_surface[2]['pdos'][0, :],
            x_axis, -data_surface[2]['pdos'][1, :],
            linewidth = 1.0,
            color = "black"
        )

        x_axis = data_surface[i]['eigen_energy']

        for orb in range(len(d_colors)):
            axs[i, 3].plot(
                x_axis,  data_surface[i]['pdos'][0, :,4 + orb],
                x_axis, -data_surface[i]['pdos'][1, :,4 + orb],
                linewidth = 1.0,
                color = d_colors[orb]
            )

        axs[i, 3].plot(
            [data_surface[i]['e_fermi'], data_surface[i]['e_fermi']],
            [np.amax(data_surface[0]['pdos'][0, :]), -np.amax(data_surface[0]['pdos'][0, :])],
            linestyle = 'dashed',
            linewidth = 1.0,
            color = 'black'
        )

        axs[i, 3].plot(
            x_axis,  data_surface[2]['pdos'][0, :],
            x_axis, -data_surface[2]['pdos'][1, :],
            linewidth = 1.0,
            color = "black"
        )

    fig_legend = np.array([["a", "b", "c", "d"], [ "e", "f", "g", "h"]])
    with np.nditer(fig_legend, flags = ['multi_index']) as iterator:
        while not iterator.finished:
            axs[iterator.multi_index[0], iterator.multi_index[1]].text(
                0.05,
                0.85,
                iterator[0],
                transform = axs[iterator.multi_index[0], iterator.multi_index[1]].transAxes
            )
            axs[iterator.multi_index[0], iterator.multi_index[1]].set(
                yticks = [],
                xlim = xlim,
                ylim = ylim
            )
            axs[iterator.multi_index[0], iterator.multi_index[1]].xaxis.set_minor_locator(AutoMinorLocator())
            iterator.iternext()

    output.savefig()
