#!/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pgf import PdfPages
import matplotlib.ticker as tck

input_files = ["nn_c5_ACF.dat", "nn_t1_ACF.dat", "nn_t5_ACF.dat"]

O_layers = [
    [1, 6, 12, 39, 44, 46, 49, 50],#   L1
    [5, 11, 43, 45],#                  L1.5
    [2, 4, 38, 40],#                   L2
    [7, 8, 9, 10, 41, 42, 47, 48],#    L2.5
    [3, 13, 26, 28],#                  L3
    [17, 18, 23, 24, 31, 32, 33, 34],# L3.5
    [14, 16, 26, 28],#                 L4
    [19, 20, 21, 22, 29, 30, 35, 36],# L4.5
    [15, 25]#                          L5
    ]

Na_layers = [
    [53, 63, 67, 68],# L1
    [52, 53, 64, 66],# L2
    [51, 57, 59, 65],# L3
    [56, 58, 60, 62],# L4
    [55, 61]#          L5
    ]

num_O_layers = len(O_layers)
num_Na_layers = len(Na_layers)

width = 0.75*(210 - 30 - 20)
height = 80

raw_charge = []
for files in input_files:
    raw_charge.append(pd.read_fwf(files, skiprows=[1, 88, 89, 90, 91]))

average_O_charge = []
average_Na_charge = []

for strain, bader_chg in enumerate(raw_charge):
    average_O_layer_charge = []
    average_Na_layer_charge = []
    for i, layer in enumerate(O_layers):
        accumulator = 0
        for j, atom in enumerate(layer):
            accumulator += bader_chg['CHARGE'][atom - 1]/len(layer)
        average_O_layer_charge.append(accumulator)
    average_O_charge.append(average_O_layer_charge)
    for i, layer in enumerate(Na_layers):
        accumulator = 0
        for j, atom in enumerate(layer):
            accumulator += bader_chg['CHARGE'][atom - 1]/len(layer)
        average_Na_layer_charge.append(accumulator)
    average_Na_charge.append(average_Na_layer_charge)

average_O_charge = np.array(average_O_charge)
c4_O_change = np.divide(average_O_charge[0], average_O_charge[1]) - 1
t5_O_change = np.divide(average_O_charge[2], average_O_charge[1]) - 1

average_Na_charge = np.array(average_Na_charge)
c4_Na_change = np.divide(average_Na_charge[0], average_Na_charge[1]) - 1
t4_Na_change = np.divide(average_Na_charge[2], average_Na_charge[1]) - 1

with PdfPages('bader_O-Na_evolution_NTOvsNNO_alt.pdf') as output:
    plt.rcParams.update({
        'font.family': 'serif',
        'font.size': 11,
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
        nrows=1,
        ncols=2,
        constrained_layout=True
        )

    axs[0].scatter(
            100*c4_Na_change,
            [1, 2, 3, 4, 5],
            marker='o',
            linestyle="dashed",
            label=r'$-5\%$',
            color='k'
            )

    axs[0].scatter(
            100*t4_Na_change,
            [1, 2, 3, 4, 5],
            marker='o',
            label=r'$5\%$',
            color='darkgray'
            )

    axs[1].scatter(
            100*c4_O_change,
            np.linspace(1, 5, 9),
            marker='o',
            label=r'$-5\%$',
            color='k'
            )

    axs[1].scatter(
            100*t5_O_change,
            np.linspace(1, 5, 9),
            marker='o',
            label=r'$5\%$',
            color='darkgray'
            )

    axs[0].set(
        ylabel='Layers',
        ylim=[6, 0],
        yticks=[],
        xlabel=r'$\Delta$Charge (ground state \%)',
        xlim=[-0.5, 0.5]
    )

    axs[1].set(
        ylim=[6, 0],
        xlim=[-0.75, 0.5],
        xlabel=r'$\Delta$Charge (ground state \%)'
    )

    axs[1].tick_params(
        axis='y',
        which='major',
        length=10.0,
        labelleft=False,
        labelright=False,
        left=False,
        right=False
    )

    axs[0].xaxis.set_minor_locator(tck.AutoMinorLocator())
    axs[1].xaxis.set_minor_locator(tck.AutoMinorLocator())

    axs[0].text(0.05, 0.90, 'a', transform=axs[0].transAxes, fontsize='large')
    axs[1].text(0.05, 0.90, 'b', transform=axs[1].transAxes, fontsize='large')
    axs[1].legend(
        ncol=2,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.4)
    )

    output.savefig()
