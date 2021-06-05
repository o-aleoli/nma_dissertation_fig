#!/usr/bin/env python3
import numpy as np
from matplotlib.backends.backend_pgf import PdfPages
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches

INPUT_FILE_NB = [
    'nno/total',
    'nno/Nb_d_sum_1',
    'nno/Nb_d_sum_2',
    'nno/Nb_d_sum_3',
    'nno/Nb_d_sum_4',
    'nno/Nb_d_sum_5'
]
INPUT_FILE_TA = [
    'nto/total',
    'nto/Ta_d_sum_1',
    'nto/Ta_d_sum_2',
    'nto/Ta_d_sum_3',
    'nto/Ta_d_sum_4',
    'nto/Ta_d_sum_5'
]
OUTPUT_FILE = 'lrpdos_tf_nno_nto'
XLIM = (-6, 9)
YLIM = (-10, 10)

DATA_NB = []
for i in INPUT_FILE_NB:
    DATA_NB.append(np.load(i+'.npz'))

DATA_TA = []
for i in INPUT_FILE_TA:
    DATA_TA.append(np.load(i+'.npz'))

MM2IN = 25.4
WIDTH = (210-30-20)/MM2IN
HEIGHT = (297-30-20)/MM2IN

"""
INPUT_DATA is a pickled, compressed file with the following numpy arrays:
    * pdos, an np.array, with dimensions (ispin, numpoints), where ispin=2 and numpoints=301 (typically)
    * eigen_energy, an np.array, with dimensions (numpoints)
    * eigen_adjust, a float
    * e_fermi, a float
"""
with PdfPages(OUTPUT_FILE+'.pdf') as OUTPUT:
    plt.rcParams.update({
        'font.family': 'serif',
        'figure.figsize': [WIDTH, HEIGHT],
        'savefig.dpi': 1000,
        'text.usetex': True,
        'pgf.rcfonts': False,
        'pgf.texsystem': 'pdflatex',
        'pgf.preamble': r'\usepackage{mathpazo,eulervm}\usepackage[utf8x]{inputenc}',
        'legend.loc': 'upper right',
        'legend.frameon': False,
        'legend.fancybox': False
    })
    FIG, AXS = plt.subplots(
        nrows=5,
        ncols=2,
        # constrained_layout=True,
    )

    x_axis = DATA_NB[0]['eigen_energy'] + DATA_NB[0]['eigen_adjust']

###############################################################################

    AXS[0, 0].plot(
        x_axis,  DATA_NB[0]['pdos'][0, ...],
        x_axis, -DATA_NB[0]['pdos'][1, ...],
        c='black',
        linewidth=1.0
    )

    AXS[0, 0].plot(
        x_axis,  DATA_NB[1]['pdos'][0, ...],
        x_axis, -DATA_NB[1]['pdos'][1, ...],
        linewidth=1.0,
        c='green'
    )

    AXS[0, 0].plot(
        [DATA_NB[0]['e_fermi'], DATA_NB[0]['e_fermi']],
        [np.amax(DATA_NB[0]['pdos'][0, ...]), -np.amax(DATA_NB[0]['pdos'][0, ...])],
        linestyle='dashed',
        linewidth=2.0,
        c='black'
    )

###############################################################################

    AXS[1, 0].plot(
        x_axis,  DATA_NB[0]['pdos'][0, ...],
        x_axis, -DATA_NB[0]['pdos'][1, ...],
        c='black',
        linewidth=1.0
    )

    AXS[1, 0].plot(
        x_axis,  DATA_NB[2]['pdos'][0, ...],
        x_axis, -DATA_NB[2]['pdos'][1, ...],
        linewidth=1.0,
        c='green'
    )

    AXS[1, 0].plot(
        [DATA_NB[0]['e_fermi'], DATA_NB[0]['e_fermi']],
        [np.amax(DATA_NB[0]['pdos'][0, ...]), -np.amax(DATA_NB[0]['pdos'][0, ...])],
        linestyle='dashed',
        linewidth=2.0,
        c='black'
    )

###############################################################################

    AXS[2, 0].plot(
        x_axis,  DATA_NB[0]['pdos'][0, ...],
        x_axis, -DATA_NB[0]['pdos'][1, ...],
        c='black',
        linewidth=1.0
    )

    AXS[2, 0].plot(
        x_axis,  DATA_NB[3]['pdos'][0, ...],
        x_axis, -DATA_NB[3]['pdos'][1, ...],
        linewidth=1.0,
        c='green'
    )

    AXS[2, 0].plot(
        [DATA_NB[0]['e_fermi'], DATA_NB[0]['e_fermi']],
        [np.amax(DATA_NB[0]['pdos'][0, ...]), -np.amax(DATA_NB[0]['pdos'][0, ...])],
        linestyle='dashed',
        linewidth=2.0,
        c='black'
    )

###############################################################################

    AXS[3, 0].plot(
        x_axis,  DATA_NB[0]['pdos'][0, ...],
        x_axis, -DATA_NB[0]['pdos'][1, ...],
        c='black',
        linewidth=1.0
    )

    AXS[3, 0].plot(
        x_axis,  DATA_NB[4]['pdos'][0, ...],
        x_axis, -DATA_NB[4]['pdos'][1, ...],
        linewidth=1.0,
        c='green'
    )

    AXS[3, 0].plot(
        [DATA_NB[0]['e_fermi'], DATA_NB[0]['e_fermi']],
        [np.amax(DATA_NB[0]['pdos'][0, ...]), -np.amax(DATA_NB[0]['pdos'][0, ...])],
        linestyle='dashed',
        linewidth=2.0,
        c='black'
    )

###############################################################################

    AXS[4, 0].plot(
        x_axis,  DATA_NB[0]['pdos'][0, ...],
        x_axis, -DATA_NB[0]['pdos'][1, ...],
        c='black',
        linewidth=1.0
    )

    AXS[4, 0].plot(
        x_axis,  DATA_NB[5]['pdos'][0, ...],
        x_axis, -DATA_NB[5]['pdos'][1, ...],
        linewidth=1.0,
        c='green'
    )

    AXS[4, 0].plot(
        [DATA_NB[0]['e_fermi'], DATA_NB[0]['e_fermi']],
        [np.amax(DATA_NB[0]['pdos'][0, ...]), -np.amax(DATA_NB[0]['pdos'][0, ...])],
        linestyle='dashed',
        linewidth=2.0,
        c='black'
    )

###############################################################################

    AXS[0, 1].plot(
        x_axis,  DATA_TA[0]['pdos'][0, ...],
        x_axis, -DATA_TA[0]['pdos'][0, ...],
        c='black',
        linewidth=1.0
    )

    AXS[0, 1].plot(
        x_axis,  DATA_TA[1]['pdos'][0, ...],
        x_axis, -DATA_TA[1]['pdos'][0, ...],
        linewidth=1.0,
        c='#B79A56'
    )

    AXS[0, 1].plot(
        [DATA_TA[0]['e_fermi'], DATA_TA[0]['e_fermi']],
        [np.amax(DATA_TA[0]['pdos'][0, ...]), -np.amax(DATA_TA[0]['pdos'][0, ...])],
        linestyle='dashed',
        linewidth=2.0,
        c='black'
    )

###############################################################################

    AXS[1, 1].plot(
        x_axis,  DATA_TA[0]['pdos'][0, ...],
        x_axis, -DATA_TA[0]['pdos'][0, ...],
        c='black',
        linewidth=1.0
    )

    AXS[1, 1].plot(
        x_axis,  DATA_TA[2]['pdos'][0, ...],
        x_axis, -DATA_TA[2]['pdos'][0, ...],
        linewidth=1.0,
        c='#B79A56'
    )

    AXS[1, 1].plot(
        [DATA_TA[0]['e_fermi'], DATA_TA[0]['e_fermi']],
        [np.amax(DATA_TA[0]['pdos'][0, ...]), -np.amax(DATA_TA[0]['pdos'][0, ...])],
        linestyle='dashed',
        linewidth=2.0,
        c='black'
    )

###############################################################################

    AXS[2, 1].plot(
        x_axis,  DATA_TA[0]['pdos'][0, ...],
        x_axis, -DATA_TA[0]['pdos'][0, ...],
        c='black',
        linewidth=1.0
    )

    AXS[2, 1].plot(
        x_axis,  DATA_TA[3]['pdos'][0, ...],
        x_axis, -DATA_TA[3]['pdos'][0, ...],
        linewidth=1.0,
        c='#B79A56'
    )

    AXS[2, 1].plot(
       [DATA_TA[0]['e_fermi'], DATA_TA[0]['e_fermi']],
        [np.amax(DATA_TA[0]['pdos'][0, ...]), -np.amax(DATA_TA[0]['pdos'][0, ...])],
        linestyle='dashed',
        linewidth=2.0,
        c='black'
    )

###############################################################################

    AXS[3, 1].plot(
        x_axis,  DATA_TA[0]['pdos'][0, ...],
        x_axis, -DATA_TA[0]['pdos'][0, ...],
        c='black',
        linewidth=1.0
    )

    AXS[3, 1].plot(
        x_axis,  DATA_TA[4]['pdos'][0, ...],
        x_axis, -DATA_TA[4]['pdos'][0, ...],
        linewidth=1.0,
        c='#B79A56'
    )

    AXS[3, 1].plot(
        [DATA_TA[0]['e_fermi'], DATA_TA[0]['e_fermi']],
        [np.amax(DATA_TA[0]['pdos'][0, ...]), -np.amax(DATA_TA[0]['pdos'][0, ...])],
        linestyle='dashed',
        linewidth=2.0,
        c='black'
    )

###############################################################################

    AXS[4, 1].plot(
        x_axis,  DATA_TA[0]['pdos'][0, ...],
        x_axis, -DATA_TA[0]['pdos'][0, ...],
        c='black',
        linewidth=1.0
    )

    AXS[4, 1].plot(
        x_axis,  DATA_TA[5]['pdos'][0, ...],
        x_axis, -DATA_TA[5]['pdos'][0, ...],
        linewidth=1.0,
        c='#B79A56'
    )

    AXS[4, 1].plot(
        [DATA_TA[0]['e_fermi'], DATA_TA[0]['e_fermi']],
        [np.amax(DATA_TA[0]['pdos'][0, ...]), -np.amax(DATA_TA[0]['pdos'][0, ...])],
        linestyle='dashed',
        linewidth=2.0,
        c='black'
    )

###############################################################################
    LINES = [
        Line2D(
            [0], [0],
            color='green',
            linewidth=1.0,
            label=r'Nb $4d$'
        ),
        Line2D(
            [0], [0],
            color='#B79A56',
            linewidth=1.0,
            label=r'Ta $5d$'
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

    FIG.legend(handles=LINES, loc='lower center', bbox_to_anchor=(0.5, 0), ncol=4)

    AXS[0, 0].set(
        xlim=XLIM,
        ylim=YLIM,
        yticks=list(YLIM)
    )

    AXS[1, 0].set(
        xlim=XLIM,
        ylim=YLIM,
        yticks=list(YLIM)
    )

    AXS[2, 0].set(
        ylabel=r'LRPDOS (eV${}^{-1}$)',
        xlim=XLIM,
        ylim=YLIM,
        yticks=list(YLIM)
    )

    AXS[3, 0].set(
        xlim=XLIM,
        ylim=YLIM,
        yticks=list(YLIM)
    )

    AXS[4, 0].set(
        xlabel='Energy (eV)',
        xlim=XLIM,
        ylim=YLIM,
        yticks=list(YLIM)
    )

    AXS[0, 1].set(
        xlim=XLIM,
        ylim=(-0.5, 0.5),
        yticks=[-0.5, 0.5]
    )

    AXS[1, 1].set(
        xlim=XLIM,
        ylim=(-0.5, 0.5),
        yticks=[-0.5, 0.5]
    )

    AXS[2, 1].set(
        xlim=XLIM,
        ylim=(-0.5, 0.5),
        yticks=[-0.5, 0.5]
    )

    AXS[3, 1].set(
        xlim=XLIM,
        ylim=(-0.5, 0.5),
        yticks=[-0.5, 0.5]
    )

    AXS[4, 1].set(
        xlabel='Energy (eV)',
        xlim=XLIM,
        ylim=(-0.5, 0.5),
        yticks=[-0.5, 0.5]
    )

    OUTPUT.savefig()