#!/usr/bin/env python3
import numpy as np
from matplotlib.backends.backend_pgf import PdfPages
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

data_Nb = [
    [np.load('nn_total_t5.npz'),  np.load('Nb_d_separated_1_t5.npz'),  np.load('Nb_d_separated_2_t5.npz')],
    [np.load('nn_total_rlx.npz'), np.load('Nb_d_separated_1_rlx.npz'), np.load('Nb_d_separated_2_rlx.npz')],
    [np.load('nn_total_c5.npz'),  np.load('Nb_d_separated_1_c5.npz'),  np.load('Nb_d_separated_2_c5.npz')],
]
d_colors = [
    '#0d8200',#dxy
    '#13bb00',#dyz
    '#ff00d9',#dz2
    '#18f500',#dxz
    '#860069'# dx2
]
legend_1 = [
    Line2D(
        [0], [0],
        color=d_colors[4],
        linewidth=1.0,
        label=r'Nb $4d_{x^2 - y^2}$'
    ),
    Line2D(
        [0], [0],
        color=d_colors[2],
        linewidth=1.0,
        label=r'Nb $4d_{3z^2 - r^2}$'
    ),
    Line2D(
        [0], [0],
        color=d_colors[0],
        linewidth=1.0,
        label=r'Nb $4d_{xy}$'
    ),
    Line2D(
        [0], [0],
        color=d_colors[1],
        linewidth=1.0,
        label=r'Nb $4d_{yz}$'
    ),
]
legend_2 = [
    Line2D(
        [0], [0],
        color=d_colors[3],
        linewidth=1.0,
        label=r'Nb $4d_{xz}$'
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
subfigure = np.array([
    ['a', 'b'],
    ['c', 'd'],
    ['e', 'f'],
])

output_file = 'nn_separated_pdos_evo'
x_lim = (1, 7)
y_lim = (-5, 5)
width = 210 - 30 - 20
height = 0.5*(297 - 30 - 20) - 15

with PdfPages(output_file + '.pdf') as output:
    plt.rcParams.update({
        'font.size': 11.0,
        'figure.figsize': [width/25.4, height/25.4],
        'savefig.dpi': 1200,
        'text.usetex': True,
        'pgf.rcfonts': False,
        'pgf.texsystem': 'pdflatex',
        'pgf.preamble': r'\usepackage{mathpazo,eulervm}\usepackage[utf8x]{inputenc}\usepackage[separate-uncertainty = true]{siunitx}',
        'legend.frameon': False,
        'legend.fancybox': False,
        'legend.fontsize': 'xx-small'
    })
    fig, axs = plt.subplots(
        nrows = 3,
        ncols = 2,
        constrained_layout = True
    )

    for strain in range(len(data_Nb)):
        x_axis = data_Nb[strain][0]['eigen_energy'] + data_Nb[strain][0]['eigen_adjust']
        for orbital in range(len(d_colors)):
            axs[strain, 0].plot(
                x_axis,  data_Nb[strain][1]['pdos'][0, :, 4 + orbital],
                x_axis, -data_Nb[strain][1]['pdos'][1, :, 4 + orbital],
                color = d_colors[orbital],
                linestyle = 'solid',
                linewidth = 1.0
            )
            axs[strain, 1].plot(
                x_axis,  data_Nb[strain][2]['pdos'][0, :, 4 + orbital],
                x_axis, -data_Nb[strain][2]['pdos'][1, :, 4 + orbital],
                color = d_colors[orbital],
                linestyle = 'solid',
                linewidth = 1.0
            )
        for layer in range(2):
            axs[strain, layer].vlines(
                data_Nb[strain][0]['e_fermi'] + data_Nb[strain][0]['eigen_adjust'],
                y_lim[0], y_lim[1],
                color = 'black',
                linestyle = 'dashed',
                linewidth = 2.0
            )
            axs[strain, layer].plot(
                x_axis,  data_Nb[strain][0]['pdos'][0, :],
                x_axis, -data_Nb[strain][0]['pdos'][1, :],
                color = 'black',
                linestyle = 'solid',
                linewidth = 1.0
            )

            axs[strain, layer].set(
                xlim = x_lim,
                ylim = y_lim,
                # xticks = range(x_lim[0], x_lim[1] + 1),
                yticks = [],
                xticklabels = []
            )

    axs[0, 0].set_title(r'$\epsilon = 5\%', x = 1, y = 1)
    axs[1, 0].set_title(r'Ground state ($\epsilon = 1\%$)', x = 1, y = 1)
    axs[2, 0].set_title(r'$\epsilon = -5\%', x = 1, y = 1)

    axs[2, 0].set(xlabel = r'Energy (\si{\electronvolt})', xticklabels = range(x_lim[0], x_lim[1] + 1))
    axs[2, 1].set(xlabel = r'Energy (\si{\electronvolt})', xticklabels = range(x_lim[0], x_lim[1] + 1))
    axs[1, 0].set(ylabel = r'pDOS (\si{\per\electronvolt})')

    axs[2, 0].legend(
        handles=legend_1,
        loc='lower center',
        bbox_to_anchor=(0.45, -1),
        ncol=3,
    )

    axs[2, 1].legend(
        handles=legend_2,
        loc='lower center',
        bbox_to_anchor=(0.3, -0.8),
        ncol=3,
    )

    iterator = np.nditer(subfigure, flags = ['multi_index'])
    for element in iterator:
        axs[iterator.multi_index].text(0.1, 0.85, element, transform = axs[iterator.multi_index].transAxes)

    output.savefig()
