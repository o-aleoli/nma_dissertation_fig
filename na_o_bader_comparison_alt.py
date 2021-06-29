#!/bin/env python3.7
import pandas as pd
import matplotlib as mpl
from matplotlib.backends.backend_pgf import PdfPages
from matplotlib import rc
import matplotlib.pyplot as plt
import matplotlib.style

NA_ION = pd.read_fwf('na_ion.dat')
O_ION = pd.read_fwf('o_ion.dat')

mm2in = 1/25.4
GEOMETRY = [0.75*(210-30-20)*mm2in, (297-30-20)/3*mm2in]

with PdfPages('na_o_bader_chg_comparison-alt.pdf') as output:
    plt.rcParams.update({
        'font.family': 'serif',
        'font.size': 11.0,
        'figure.figsize': GEOMETRY,
        'savefig.dpi': 96,
        'text.usetex': True,
        'pgf.rcfonts': False,
        'pgf.texsystem': 'pdflatex',
        'pgf.preamble': r'\usepackage{mathpazo,eulervm}',
        'legend.frameon': False,
        'legend.fancybox': False
        })

    fig, (axs0, axs1) = plt.subplots(
            nrows=1,
            ncols=2,
            sharex=True,
            constrained_layout=True
    )

    axs0.scatter(
            NA_ION['NNO'],
            NA_ION['layer'],
            c='black',
            marker='o',
            label='NNO'
    )

    axs1.scatter(
            O_ION['NNO'],
            O_ION['layer'],
            c='black',
            marker='o',
            label='NNO'
    )

    axs0.set(
        ylabel='Layers',
        ylim=[5.5, 0.5],
        yticks=[],
        xlabel=r'$\Delta$Charge ($\%$ bulk)'
    )

    axs1.set(
        ylim=[5.5, 0.5],
        xlabel=r'$\Delta$Charge ($\%$ bulk)'
    )

    axs1.tick_params(
        axis='y',
        which='major',
        length=10.0,
        labelleft=False,
        labelright=False,
        left=False,
        right=False
    )

    axs0.text(0.05, 1.025, '(a)', transform=axs0.transAxes, fontsize='large')
    axs1.text(0.05, 1.025, '(b)', transform=axs1.transAxes, fontsize='large')

    output.savefig()
