#!/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pgf import PdfPages
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

nn_layers = [
    [69, 71, 85, 86],  # L1
    [70, 72, 82, 84],  # L2
    [73, 75, 81, 83],  # L3
    [74, 76, 78, 80],  # L4
    [77, 79]  # L5
]
nt_layers = [
    [73, 77, 82, 86],
    [69, 72, 78, 81],
    [74, 76, 83, 85],
    [70, 71, 79, 80],
    [75, 84]
]
nn_bulk = range(17, 20 + 1)
nt_bulk = range(65, 80 + 1)
geometry = [0.5 * (210 - 30 - 20) / 25.4, 80 / 25.4]

nn_bader = pd.read_fwf('nn_rlx_ACF.dat', skiprows=[1, 88, 89, 90, 91])
nt_bader = pd.read_fwf('nt_rlx_ACF.dat', skiprows=[1, 88, 89, 90, 91])
nn_bulk_bader = pd.read_fwf('nn_bulk_ACF.dat', skiprows=[1, 22, 23, 24, 25])
nt_bulk_bader = pd.read_fwf('nt_bulk_ACF.dat', skiprows=[1, 82, 83, 84, 85])

# averages individual charges for each layer
nn_layer_average = []
for _, layer in enumerate(nn_layers):
    accumulator = 0
    for _, atom in enumerate(layer):
        accumulator += nn_bader['CHARGE'][atom - 1]/len(layer)
    nn_layer_average.append(accumulator)
nn_average_charge = np.array(nn_layer_average)

nn_bulk_average_charge = 0
for atom in nn_bulk:
    nn_bulk_average_charge += nn_bulk_bader['CHARGE'][atom - 1]/len(nn_bulk)

nt_layer_average = []
for _, layer in enumerate(nt_layers):
    accumulator = 0
    for _, atom in enumerate(layer):
        accumulator += nt_bader['CHARGE'][atom - 1]/len(layer)
    nt_layer_average.append(accumulator)
nt_average_charge = np.array(nn_layer_average)

nt_bulk_average_charge = 0
for atom in nt_bulk:
    nt_bulk_average_charge += nt_bulk_bader['CHARGE'][atom - 1]/len(nt_bulk)

nn_bader_charge_change = 100*(nn_average_charge/nn_bulk_average_charge - 1)
nt_bader_charge_change = 100*(nt_average_charge/nt_bulk_average_charge - 1)

with PdfPages('nn_nt_bader.pdf') as output:
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
        nn_average_charge,
        [1, 2, 3, 4, 5],
        marker='s',
        label=r'NaNbO$_3$',
        color='#4CB276'
    )
    axs.scatter(
        nt_average_charge,
        [1, 2, 3, 4, 5],
        marker='x',
        label=r'NaTaO$_3$',
        color="#B79A56"
    )

    axs.legend(
        loc = 'lower center',
        bbox_to_anchor = (0.5, -0.35),
        ncol=2
    )

    axs.set(
        xlabel=r'Total charge (electrons)',
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

    axs.xaxis.set_minor_locator(AutoMinorLocator())

    output.savefig()
