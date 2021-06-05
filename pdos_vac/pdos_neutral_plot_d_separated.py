#!/usr/bin/env python3
import numpy as np
from matplotlib.backends.backend_pgf import PdfPages
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

input_files_Na0 = [
    'Na0_total',
    'Na0_Nb_d_separated',
    'Na0_O_p_sum'
]
input_files_O0 = [
   'O0_total',
   'O0_Nb_d_separated',
   'O0_O_p_sum'
]
output_file = 'pdos_neutral_d_separated_333'
xlim = (-4.25, 9)
ylim = (-400, 400)

data_Na0 = []
for i in input_files_Na0:
    data_Na0.append(np.load(i+'.npz'))

data_O0 = []
for i in input_files_O0:
   data_O0.append(np.load(i+'.npz'))

width = 210 - 30 - 20
height = 80

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
        'figure.figsize': [width/25.4, height/25.4],
        'savefig.dpi': 1200,
        'text.usetex': True,
        'pgf.rcfonts': False,
        'pgf.texsystem': 'pdflatex',
        'pgf.preamble': r'\usepackage{mathpazo,eulervm}\usepackage[utf8x]{inputenc}',
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

###############################################################################

    x_axis = data_O0[0]['eigen_energy']

    ax_inset = zoomed_inset_axes(
        axs[0],
        3.5,
        loc = 'center',
        bbox_to_anchor = [0.475, 0.175],
        bbox_transform = axs[0].transAxes,
        axes_kwargs = {'box_aspect': 1.0}
    )

    axs[0].plot(
        x_axis,  data_O0[0]['pdos'][0, ...],
        x_axis, -data_O0[0]['pdos'][0, ...],
        color = 'black',
        linewidth = 1.0
    )

    for orb in range(len(d_colors)):
        axs[0].plot(
            x_axis,  data_O0[1]['pdos'][0, ...,4 + orb],
            x_axis, -data_O0[1]['pdos'][0, ...,4 + orb],
            linewidth = 1.0,
            color = d_colors[orb]
        )

    axs[0].plot(
        x_axis,  data_O0[2]['pdos'][0, ...],
        x_axis, -data_O0[2]['pdos'][0, ...],
        linewidth = 1.0,
        color = 'red'
    )

    axs[0].plot(
        [data_O0[0]['e_fermi'], data_O0[0]['e_fermi']],
        [np.amax(data_O0[0]['pdos'][0, ...]), -np.amax(data_O0[0]['pdos'][0, ...])],
        linestyle = 'dashed',
        color = 'black',
        linewidth = 2.0
    )

    ax_inset.plot(
        x_axis,  data_O0[0]['pdos'][0, ...],
        x_axis, -data_O0[0]['pdos'][0, ...],
        color = 'black',
        linewidth = 1.0
    )

    for orb in range(len(d_colors)):
        ax_inset.plot(
            x_axis,  data_O0[1]['pdos'][0, ...,4 + orb],
            x_axis, -data_O0[1]['pdos'][0, ...,4 + orb],
            linewidth = 1.0,
            color = d_colors[orb]
        )

    ax_inset.plot(
        x_axis,  data_O0[2]['pdos'][0, ...],
        x_axis, -data_O0[2]['pdos'][0, ...],
        linewidth = 1.0,
        color = 'red'
    )

    ax_inset.plot(
        [data_O0[0]['e_fermi'], data_O0[0]['e_fermi']],
        [np.amax(data_O0[0]['pdos'][0, ...]), -np.amax(data_O0[0]['pdos'][0, ...])],
        linestyle = 'dashed',
        linewidth = 1.0,
        color = 'black'
    )

###############################################################################

    x_axis = data_Na0[0]['eigen_energy']

    axs[1].plot(
        x_axis,  data_Na0[0]['pdos'][0, ...],
        x_axis, -data_Na0[0]['pdos'][0, ...],
        c='black',
        linewidth=1.0
    )

    for orb in range(len(d_colors)):
        axs[1].plot(
            x_axis,  data_Na0[1]['pdos'][0, ...,4 + orb],
            x_axis, -data_Na0[1]['pdos'][0, ...,4 + orb],
            linewidth=1.0,
            c=d_colors[orb]
        )

    axs[1].plot(
        x_axis,  data_Na0[2]['pdos'][0, ...],
        x_axis, -data_Na0[2]['pdos'][0, ...],
        linewidth=1.0,
        c='red'
    )

    axs[1].plot(
        [data_Na0[0]['e_fermi'], data_Na0[0]['e_fermi']],
        [np.amax(data_Na0[0]['pdos'][0, ...]), -np.amax(data_Na0[0]['pdos'][0, ...])],
        linestyle='dashed',
        linewidth=2.0,
        c='black'
    )

###############################################################################

    LINES = [
        Line2D(
            [0], [0],
            color="#860069",
            linewidth=1.0,
            label=r'Nb $4d_{x^2 - y^2}$'
        ),
        Line2D(
            [0], [0],
            color="#ff00d9",
            linewidth=1.0,
            label=r'Nb $4d_{3z^2 - r^2}$'
        ),
        Line2D(
            [0], [0],
            color="#0d8200",
            linewidth=1.0,
            label=r'Nb $4d_{xy}$'
        ),
        Line2D(
            [0], [0],
            color="#18f500",
            linewidth=1.0,
            label=r'Nb $4d_{xz}$'
        ),
        Line2D(
            [0], [0],
            color="#13bb00",
            linewidth=1.0,
            label=r'Nb $4d_{zy}$'
        ),
        Line2D(
            [0], [0],
            color='red',
            linewidth=1.0,
            label=r'O $2p$'
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
        ),
    ]

    axs[1].legend(
        handles = LINES,
        loc = 'lower center',
        bbox_to_anchor = (0.5, -1),
        ncol = 5
    )

    axs[0].set(
        xlim = xlim,
        ylim = ylim,
        yticks = [],
        xticklabels = [],
        ylabel = r'pDOS (eV${}^{-1}$)'
    )

    axs[1].set(
        xlim = xlim,
        ylim = ylim,
        yticks = [],
        xlabel = 'Energy (eV)',
        ylabel = r'pDOS (eV${}^{-1}$)'
    )

    ax_inset.set(
        xlim = (3.1, 3.7),
        ylim = (-70, 70),
        xticks = [],
        yticks = []
    )

    mark_inset(
        axs[0],
        ax_inset,
        loc1 = 2,
        loc2 = 4,
        edgecolor = 'gray'
    )

    axs[0].text(
        0.01,
        0.85,
        "a",
        transform = axs[0].transAxes
    )

    axs[1].text(
        0.01,
        0.85,
        "b",
        transform = axs[1].transAxes
    )

    output.savefig()
