#!/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pgf import PdfPages

input_files = ['nn_c5_ACF.dat', 'nn_t1_ACF.dat', 'nn_t5_ACF.dat']
#BADER_NTO = []
nn_layers = [
    [69, 71, 85, 86],  # L1
    [70, 72, 82, 84],  # L2
    [73, 75, 81, 83],  # L3
    [74, 76, 78, 80],  # L4
    [77, 79]  # L5
]
# LAYERS_NTO = [
#    [73, 77, 82, 86],
#    [69, 72, 78, 81],
#    [74, 76, 83, 85],
#    [70, 71, 79, 80],
#    [75, 84]
#    ]
num_layers_nn = len(nn_layers)
#NUMLAYERS_NTO = len(LAYERS_NTO)
geometry = [0.5*(210-30-20)/25.4, 80/25.4]

nn_bader = []
for files in input_files:
    nn_bader.append(pd.read_fwf(files, skiprows=[1, 88, 89, 90, 91]))

#BADER_BULK_NNO = pd.read_fwf('../../../bulk/1x1x1/bader/ACF.dat', skiprows=[1, 22, 23, 24, 25])

#AVERAGE_BULK_NNO = np.average(np.array(BADER_BULK_NNO))
nn_average_charge = []

# averages individual charges for each layer
for strain, bader_chg in enumerate(nn_bader):
    layer_average = []
    for i, layer in enumerate(nn_layers):
        accumulator = 0
        for j, atom in enumerate(layer):
            accumulator += bader_chg['CHARGE'][atom - 1]/len(layer)
        layer_average.append(accumulator)
    nn_average_charge.append(layer_average)

nn_average_valence = np.array(nn_average_charge)
nn_c4_change = np.divide(nn_average_charge[0], nn_average_charge[1]) - 1
nn_t4_change = np.divide(nn_average_charge[2], nn_average_charge[1]) - 1
#BADER_CHANGE_NNO = nn_average_valence/BADER_BULK_NNO['CHARGE'][19] - 1

# for files in input_files:
#    BADER_NTO.append(pd.read_fwf('../../../../o-NaTaO3/bandGapEngineering/' + files + '/baderCharge/ACF.dat', skiprows=[1, 88, 89, 90, 91]))
#
#BADER_BULK_NTO = 8.175205
#
#AVERAGE_CHARGE_NTO = []
#
# for strain, bader_chg in enumerate(BADER_NTO):
#    layer_average = []
#    for i, layer in enumerate(LAYERS_NTO):
#        accumulator = 0
#        for j, atom in enumerate(layer):
#            accumulator += bader_chg['CHARGE'][atom - 1]/len(layer)
#        layer_average.append(accumulator)
#    AVERAGE_CHARGE_NTO.append(layer_average)
#
#AVERAGE_CHARGE_NTO = np.array(AVERAGE_CHARGE_NTO)
#BADER_CHANGE_NTO = AVERAGE_CHARGE_NTO/BADER_BULK_NTO - 1


with PdfPages('bader_strain_evolution_NTOvsNNO-alt.pdf') as output:
    plt.rcParams.update({
        'font.family': 'serif',
        'figure.figsize': geometry,
        'savefig.dpi': 1200,
        'text.usetex': True,
        'pgf.rcfonts': False,
        'pgf.texsystem': 'pdflatex',
        'pgf.preamble': r'\usepackage{mathpazo,eulervm}\usepackage[utf8x]{inputenc}',
        'legend.frameon': False,
        'legend.fancybox': False,
        'font.size': 11
    })

    fig, axs = plt.subplots(constrained_layout=True)

    axs.scatter(
        100*nn_c4_change,
        [1, 2, 3, 4, 5],
        marker='o',
        label=r'$-5\%$',
        color='k'
    )
    axs.scatter(
        100*nn_t4_change,
        [1, 2, 3, 4, 5],
        marker='o',
        label=r'$5\%$',
        color='darkgray'
    )
#    axs.plot(
#            [1, 2, 3, 4, 5],
#            100*BADER_CHANGE_NNO[2, ...],
#            marker='o',
#            label=r'$5\%$',
#            color='lightgray'
#            )
#
#    axs.plot(
#            [1, 2, 3, 4, 5],
#            100*BADER_CHANGE_NTO[0, ...],
#            marker='o',
#            label=r'$-5\%$',
#            color='k'
#            )
#    axs[1].plot(
#            [1, 2, 3, 4, 5],
#            100*BADER_CHANGE_NTO[1, ...],
#            marker='o',
#            label=r'$0\%$',
#            color='darkgray'
#            )
#    axs[1].plot(
#            [1, 2, 3, 4, 5],
#            100*BADER_CHANGE_NTO[2, ...],
#            marker='o',
#            label=r'$5\%$',
#            color='lightgray'
#            )

    axs.legend(loc="lower left")

    axs.set(
        xlabel=r'$\Delta$Charge (ground state \%)',
        xlim=[-2.0, 2.0],
        ylabel=r'Layers',
        ylim=[5.5, 0.5]
    )

    axs.tick_params(
        axis='y',
        which='major',
        length=10.0,
        labelleft=False,
        labelright=False,
        left=False,
        right=False
    )
#    axs.set_xticks([1, 2, 3, 4, 5])

#    axs[0].annotate(r'NNO', xy=(0.35, 1.05), xycoords='axes fraction')
#    axs[1].annotate(r'NTO', xy=(0.35, 1.05), xycoords='axes fraction')

    output.savefig()
